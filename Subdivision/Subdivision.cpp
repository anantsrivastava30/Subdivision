#include <fstream>
#include <vector>
#include "OBJFileReader.h"
#include "Solid.h"
#include "iterators.h"
#include "SolidDelegate.h"
#include <iostream>
#include <stdio.h>
#include <limits.h>
#include "mst.h"

using namespace std;
using namespace MeshLib;
using namespace mst;


int main(int argc, char *argv[])
{
    cout<< "HEYYYYY" << endl;
    // Read in the obj file
    Solid mesh;
    OBJFileReader of;
    std::ifstream in(argv[1]);
    of.readToSolid(&mesh, in);
    
    /******************* Put you FD processing here *********************/
    
    cout<< "Hello"<<endl;
    Solid newMesh;
    SolidDelegate delegate;
    
    cout<< "New Mesh"<<endl;
    mesh.UpdateNormals();
    //Dual Mesh
    ofstream fp;
    fp.open ("fid_centroid.txt");
    
    //dual vertices
    int N = mesh.numFaces();
    vector< vector<double> > Centroid;
    
    int k = 0;
    for (SolidFaceIterator fiter(&mesh); !fiter.end(); ++fiter)
    {
        Face *f = *fiter;
        HalfEdge *he = f->halfedge();
        Vertex *vertices[3];
        vertices[0] = he->source();
        vertices[1] = he->he_next()->source();
        vertices[2] = he->he_next()->he_next()->source();
        
        Vertex *v = delegate.createVertex(&newMesh, newMesh.numVertices() + 1);
        v->point() = (vertices[0]->point() + vertices[1]->point() + vertices[2]->point())/3.0;
        
        Centroid.push_back(vector<double>());
        Centroid[k].push_back(f->id());
        Centroid[k].push_back(v->point()[0]);
        Centroid[k].push_back(v->point()[1]);
        Centroid[k].push_back(v->point()[2]);
        fp << f->id() << "\t" << v->point()[0] << "\t" << v->point()[1] << "\t" << v->point()[2] <<"\n";
        cout<< "Vertex Found!"<< Centroid[k][1]<< Centroid[k][2] << Centroid[k][3] << endl;
        k++;
    }
    
    fp.close();
    //dual faces
    k = 0;
    vector< vector<double> > FaceVertex;
    ofstream fout("cutting_edge.txt");
    int fid;
    for (SolidFaceIterator fiter(&mesh); !fiter.end(); ++fiter)
    {
        
        Face *f = *fiter;
        HalfEdge * he = f->halfedge();
        
        for(int i=0;i<3;i++)
        {
            Vertex *v1 = he->source();
            Vertex *v2 = he->he_sym()->source();
            fid = he->he_sym()->face()->id();
            fout << f->id() << "  " << fid << "  " << v1->id() << "  " << v2->id()  << endl;
            FaceVertex.push_back(vector<double>());
            FaceVertex[k].push_back(f->id());
            FaceVertex[k].push_back(fid);
            FaceVertex[k].push_back(v1->id());
            FaceVertex[k].push_back(v2->id());
            cout <<	FaceVertex[k][0] << " adj face : " << FaceVertex[k][1] << endl;
            he = he->he_next();
            k++;
        }
    }
    fout.close();
    
    //make adjancy matrix for dual
    vector< vector<double> > V(N, vector<double>(N, 0));
    for (int i = 0; i < FaceVertex.size(); i++)
    {
        int v1;
        int v2;
        double dist;
        v1 = FaceVertex[i][0];
        v2 = FaceVertex[i][1];
        if (Centroid[v1 - 1][0] != v1 && Centroid[v2 - 1][0] != v2)
        {
            cout << "discrepancy" << endl;
        }
        dist = sqrt(pow((Centroid[v1 - 1][1] - Centroid[v2 - 1][1]), 2) + pow((Centroid[v1 - 1][2] - Centroid[v2 - 1][2]), 2) + pow((Centroid[v1 - 1][3] - Centroid[v2 - 1][3]), 2));
        V[v1 - 1][v2 - 1] = dist;
        V[v2 - 1][v1 - 1] = dist;
    }
    cout << V.size() << V[0].size()<< endl;
    //MST
    vector < vector < double > > X;
    X.push_back(vector < double >());
    X[0].push_back(0);
    X[0].push_back(2);
    X[0].push_back(0);
    X[0].push_back(6);
    X[0].push_back(0);
    X.push_back(vector < double >());
    X[1].push_back(2);
    X[1].push_back(0);
    X[1].push_back(3);
    X[1].push_back(8);
    X[1].push_back(5);
    X.push_back(vector < double >());
    X[2].push_back(0);
    X[2].push_back(3);
    X[2].push_back(0);
    X[2].push_back(0);
    X[2].push_back(7);
    X.push_back(vector < double >());
    X[3].push_back(6);
    X[3].push_back(8);
    X[3].push_back(0);
    X[3].push_back(0);
    X[3].push_back(9);
    X.push_back(vector < double >());
    X[4].push_back(0);
    X[4].push_back(5);
    X[4].push_back(7);
    X[4].push_back(9);
    X[4].push_back(0);
    
    //{ {0, 2, 0, 6, 0},
    //{ 2, 0, 3, 8, 5 },
    //{ 0, 3, 0, 0, 7 },
    //{ 6, 8, 0, 0, 9 },
    //{ 0, 5, 7, 9, 0 },
    //};
    MST tree;
    vector < int > prim = tree.primMST(V);
    cout << prim[0] << endl;
    
    vector < vector <double > > corel;
    for (int i = 0; i < prim.size(); i++)
    {
        corel.push_back(vector < double >());
        corel[i].push_back(i+1);
        for (int k = 0; k < 3; k++)
        {
            corel[i].push_back(FaceVertex[i*3 + k][1]);
        }
    }
    
    // matrix representation
    vector < vector <int > > cm(mesh.numVertices() + 1, vector < int>(mesh.numVertices() + 1, 0));
    
    for (int i = 0; i < prim.size(); i++)
    {
        if (prim[i] == -1)
        {
            cout << "originating vertex : " << i+1 << endl;
        }
        if (prim[i] + 1 != corel[i][1])
        {
            cm[FaceVertex[(i) * 3 + 0][2]][FaceVertex[(i) * 3 + 0][3]] = 1; cm[FaceVertex[(i) * 3 + 0][3]][FaceVertex[(i) * 3 + 0][2]] = 1;
            cout << FaceVertex[(i)*3 + 0][2] << " " << FaceVertex[(i) * 3 + 0][3] << endl;
        }
        if (prim[i] + 1 != corel[i][2])
        {
            cm[FaceVertex[(i) * 3 + 1][2]][FaceVertex[(i) * 3 + 1][3]] = 1;	cm[FaceVertex[(i) * 3 + 1][3]][FaceVertex[(i) * 3 + 1][2]] = 1;
            cout << FaceVertex[(i) * 3 + 1][2] << " " << FaceVertex[(i) * 3 + 1][3] << endl;
        }
        if (prim[i] + 1 != corel[i][3])
        {
            cm[FaceVertex[(i) * 3 + 0][2]][FaceVertex[(i) * 3 + 0][3]] = 1; cm[FaceVertex[(i) * 3 + 0][3]][FaceVertex[(i) * 3 + 0][2]] = 1;
            cout << FaceVertex[(i) * 3 + 2][2] << " " << FaceVertex[(i) * 3 + 2][3] << endl;
        }
    }
    //cout << " checker " <<cm[634][2418] << endl;
    
    
    
    // matrix representation pruning
    int probe = 1;
    int count = 0;
    int x, y;
    while (probe != 0)
    {
        probe = 0;
        for (int i=0; i < cm.size(); i++)
        {
            count = 0;
            for (int j=0; j < cm[i].size(); j++)
            {
                if (cm[i][j] == 1)
                {
                    x = i; y = j;
                    count++;
                }
            }
            /*if (count > 1)
             {
             probe = 0;
             }*/
            if (count == 1)
            {
                cm[x][y] = 0; cm[y][x] = 0;
                probe ++;
            }
            //cout << "sequence :" << probe << endl;
        }
        cout << "hallo probe ! " << probe << endl;
    }
    
    
    
    
    
    
    
    
    
    // Write out the resultant obj file
    int vObjID = 1;
    std::map<int, int> vidToObjID;
    
    std::ofstream os(argv[2]);
    
    SolidVertexIterator iter(&mesh);
    
    for(; !iter.end(); ++iter)
    {
        Vertex *v = *iter;
        Point p = v->point();
        os << "Vertex " << v->id() << " " << p[0] << " " << p[1] << " " << p[2] << std::endl;
        vidToObjID[v->id()] = vObjID++;
    }
    os << "# " << (unsigned int)mesh.numVertices() << " vertices" << std::endl;
    
    /*float u = 0.0, v = 0.0;
     for(iter.reset(); !iter.end(); ++iter)
     {
     Vertex *vv = *iter;
     std::string key( "uv" );
     std::string s = Trait::getTraitValue (vv->string(), key );
     if( s.length() > 0 )
     {
     sscanf( s.c_str (), "%f %f", &u, &v );
     }
     os << "vt " << u << " " << v << std::endl;
     }
     os << "# " << (unsigned int)mesh.numVertices() << " texture coordinates" << std::endl;
     */
    SolidFaceIterator fiter(&mesh);
    for(; !fiter.end(); ++fiter)
    {
        Face *f = *fiter;
        FaceVertexIterator viter(f);
        os << "Face " << f->id() << " " ;
        for(; !viter.end(); ++viter)
        {
            Vertex *v = *viter;
            os << vidToObjID[v->id()] << " ";
        }
        os << std::endl;
    }
    
    /*for (int i = 1; i < prim.size(); i++)
     {
     if (prim[i]+1 == corel[i][0])
     {
     os << "Edge " << FaceVertex[(i - 1) * 3 + 0][2] << " " << FaceVertex[(i - 1) * 3 + 0][3] << " {sharp}"<< endl;
     }
     if (prim[i]+1 == corel[i][1])
     {
     os << "Edge " << FaceVertex[(i - 1) * 3 + 1][2] << " " << FaceVertex[(i - 1) * 3 + 1][3] << " {sharp}" << endl;
     }
     if (prim[i]+1 == corel[i][2])
     {
     os << "Edge " << FaceVertex[(i - 1) * 3 + 2][2] << " " << FaceVertex[(i - 1) * 3 + 2][3] << " {sharp}" << endl;
     }
     }*/
    count = 0;
    for (int i = 0; i < cm.size(); i++)
    {
        for (int j = 0; j < cm[i].size(); j++)
        {
            if (cm[i][j] == 1)
            {   count++;
                os << "Edge " << i << " " << j << " {sharp}" << endl;
            }
        }
    }
    cout<< count<<endl;
    
    os.close();
    
    return 0;
}