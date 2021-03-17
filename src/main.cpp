//
//  main.cpp
//  3DA_project_CageBasedDef_bin
//
//  Created by Benjamin Barral on 05/02/2019.
//

#include <stdio.h>
#include <cstdlib>
#define DIM_X 0
#define DIM_Y 1
#define DIM_Z 2

#include <iostream>
#include <igl/readPLY.h>
#include <igl/writePLY.h>
#include <igl/writeOBJ.h>
#include <igl/file_exists.h>
#include <Eigen/Geometry>
#include <cmath>
#include <iterator>
#include "MeanValueCoordController.hpp"
#include "MeshProcessor.hpp"
#include <random>
#include "CageGenerator.hpp"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/avg_edge_length.h>

using namespace std;
using namespace Eigen;
using namespace igl::opengl::glfw;

string mesh_file_name = "../../Horse_2.obj";

float sparseness_cage = 0.5; // For automatic cage generation : CHANGE THIS for a sparser or denser cage : the larger the parameter the denser the cage

void cage_deform1(MatrixXd & V_cage, const MatrixXi & F_cage, Vector3d cage_barycenter){
    int num_vertices_cage = V_cage.rows();
    double avg_edge_len = igl::avg_edge_length(V_cage, F_cage);
    Vector3d dir = Vector3d(0, 0, 0);
    double res = 0;
    for(int i = 0; i < 15; i++){
        double theta = acos((double)rand() / RAND_MAX * 2.0 - 1.0);
        double phi = 2.0 * ((double)rand() / RAND_MAX);
        Vector3d dir_cand = Vector3d(sin(theta) * cos(phi), cos(theta), sin(theta) * cos(phi));
        double value = 0;
        for(int j = 0; j < num_vertices_cage; j++){
            Vector3d c2v = Vector3d(V_cage.row(j)) - cage_barycenter;
            double cutoff = c2v.normalized().dot(dir_cand);
            if (cutoff > 0.95){
                value += c2v.dot(dir_cand);
            }
        }
        if (value > res){
            res = value;
            dir = dir_cand;
        }
    }

    double theta = acos(1.0 -  0.01 * (double)rand() / RAND_MAX);
    double phi = 2.0 * ((double)rand() / RAND_MAX);
    Vector3d dir_cand = Vector3d(sin(theta) * cos(phi), cos(theta), sin(theta) * cos(phi));
    Matrix3d R;
    R = Quaterniond().setFromTwoVectors(Vector3d(0,1,0),dir_cand);

    for(int i = 0; i < num_vertices_cage; i++){
        Vector3d c2v = Vector3d(V_cage.row(i)) - cage_barycenter;
        double cutoff = c2v.normalized().dot(dir);
        if (cutoff > 0.95){
            V_cage.row(i) = R * c2v + cage_barycenter;
        }
    }
}

void cage_deform2(MatrixXd & V_cage, const MatrixXi & F_cage, Vector3d cage_barycenter){
    int num_vertices_cage = V_cage.rows();
    double avg_edge_len = igl::avg_edge_length(V_cage, F_cage);
    int seed = rand() % num_vertices_cage;
    double theta = acos(1.0 -  0.01 * (double)rand() / RAND_MAX);
    double phi = 2.0 * ((double)rand() / RAND_MAX);
    Vector3d dir_cand = Vector3d(sin(theta) * cos(phi), cos(theta), sin(theta) * cos(phi));
    Matrix3d R;
    R = Quaterniond().setFromTwoVectors(Vector3d(0,1,0),dir_cand);
    for(int i = 0; i < num_vertices_cage; i++){
        Vector3d c2v = Vector3d(V_cage.row(i)) - cage_barycenter;
        double d = Vector3d(V_cage.row(i) - V_cage.row(seed)).norm();
        if (d < 2 * avg_edge_len){
            V_cage.row(i) = R * c2v + cage_barycenter;
        }
    }
}

int main(int argc, char *argv[])
{
    srand(time(NULL));
    // LOAD MESHES
    MatrixXd V_mesh,V_cage;
    MatrixXi F_mesh,F_cage;
    if (mesh_file_name.find("obj") != string::npos) igl::readOBJ(mesh_file_name,V_mesh,F_mesh);
    else if (mesh_file_name.find("ply") != string::npos) igl::readPLY(mesh_file_name,V_mesh,F_mesh);
    else if (mesh_file_name.find("off") != string::npos) igl::readOFF(mesh_file_name,V_mesh,F_mesh);
    else {
        cout << "Mesh file not recognized " << endl;
        return 0;
    }
    
    // Generate automatic cage
    MatrixXd V_cage_automatic, V_cage_automatic_smooth;
    float lambda_smooth_implicit = .6;
    int num_iterations_smoothing = 2;
    CageGenerator cage_generator(V_mesh, F_mesh, num_iterations_smoothing, lambda_smooth_implicit, sparseness_cage);

    // Generate cage
    MatrixXd bb_vertices;
    MatrixXi bb_faces;
    vector<MatrixXd> feat_voxels_vertices;
    vector<MatrixXi> feat_voxels_faces;
    cage_generator.ComputeCage(bb_vertices, bb_faces, feat_voxels_vertices, feat_voxels_faces);
    
    //V_cage_automatic = cage_generator.GetCage();
    V_cage_automatic_smooth = cage_generator.GetSmoothCage();
    MatrixXi F_cage_automatic = cage_generator.GetCageFaces();
    
    // Show automatically generated cage
    V_cage = V_cage_automatic_smooth;
    F_cage = F_cage_automatic;
    
    // Get barycenter and extreme points of mesh
    MeshProcessor cage_processor(V_cage,F_cage);
    Vector3d cage_barycenter = cage_processor.GetBarycenter();
    MeshProcessor mesh_processor(V_mesh,F_mesh);
    Vector3d mesh_barycenter = mesh_processor.GetBarycenter();
    int dim = DIM_X;
    int max_vertY_ind = cage_processor.GetMaximumVertexIndex(dim);
    double bb_size = mesh_processor.GetBoundingBoxSize();
    
    // Print out number of triangles
    int num_faces_mesh = F_mesh.rows();
    int num_vertices_cage = V_cage.rows();
    int num_faces_cage = F_cage.rows();
    int num_vertices_mesh = V_mesh.rows();
    cout << "Cage : " << num_faces_cage << " triangles, " << num_vertices_cage << " vertices." << endl;
    cout << "Mesh : " << num_faces_mesh << " triangles, " << num_vertices_mesh << " vertices." << endl;
            
    // Mean Value Coordinates
    MeanValueCoordController mVCoord_controller(V_mesh,V_cage,F_mesh,F_cage, bb_size);
    mVCoord_controller.ComputeMVWeights();
    MatrixXd V_mesh_deformed = V_mesh;
    V_mesh_deformed = mVCoord_controller.MVInterpolate();
    MatrixXd V_cage_deformed = V_cage;
    mVCoord_controller.SetDeformedCage(V_cage_deformed);

    cage_deform1(V_cage_deformed, F_cage, cage_barycenter);
    //cage_deform2(V_cage, F_cage, cage_barycenter);

    mVCoord_controller.SetDeformedCage(V_cage_deformed);
    V_mesh_deformed = mVCoord_controller.MVInterpolate();

    igl::writeOBJ("../../deform1.obj", V_mesh_deformed, F_mesh);    
}
