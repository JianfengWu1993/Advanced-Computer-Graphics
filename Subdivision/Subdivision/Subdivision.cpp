#include <fstream>
#include <vector>

#include "OBJFileReader.h"
#include "Solid.h"
#include "iterators.h"
#include "SolidDelegate.h"
#include <math.h>
#define PI 3.14159265
using namespace MeshLib;



void main(int argc, char *argv[])
{
	// Read in the obj file
	Solid mesh;
	OBJFileReader of;
	std::ifstream in(argv[1]);
	of.readToSolid(&mesh, in);

	/******************* Put you subdivision processing here *********************/

	Solid nmesh;
	int numv = mesh.numVertices();
	int numf = 1;
	for (SolidFaceIterator faceiter(&mesh); !faceiter.end(); ++faceiter) 
	{
		SolidDelegate d;
		Face * f = * faceiter;
		HalfEdge * he = f->halfedge();
		Vertex * v[3];

		for (int i = 0; i < 3; i++)
		{
			v[i] = he->vertex();
			he = he->he_next();
			if (!nmesh.idVertex(v[i]->id()))
			{
				Vertex * nv = d.createVertex(&nmesh,v[i]->id());

				if(!v[i]->boundary())
				{
					Point sum;
					int n = 0;
					for (VertexEdgeIterator veiter(v[i]); !veiter.end(); ++veiter)
					{
						Vertex * vi;
						Solid::tEdge e = *veiter;
						HalfEdge *hf = e->halfedge(0);

						if (hf->target() == v[i])
							vi = hf->source();
						else
							vi = hf->target();
						sum += vi->point();
						n++;
					}
					double alpha = (0.625 - pow((3 / 8 + cos(2 * PI / n) * 1 / 4), 2)) / n;

					nv->point() = v[i]->point()*(1 - n*alpha) + sum*alpha;
				}
				else
				{
					Point pf, pb;
					for (VertexInHalfedgeIterator vihf(&mesh,v[i]); !vihf.end(); ++vihf)
					{
						HalfEdge *hf = *vihf;
						if (hf->source()->boundary())
						{
							pf = hf->source()->point();
							break;
						}
					}
					for (VertexOutHalfedgeIterator vohf(&mesh,v[i]); !vohf.end(); ++vohf)
					{
						HalfEdge *hf = *vohf;
						if (hf->target()->boundary())
						{
							pb = hf->target()->point();
							break;
						}
					}
					nv->point() = v[i]->point()*0.75 + pf*0.125 + pb*0.125;
				}

			}		

		}
		
		
		Vertex *nv01, *nv12, *nv20;
		Point pt, pb;
		Edge * e = mesh.vertexEdge(v[0], v[1]);
		if (e->string().empty())
		{
			nv01 = d.createVertex(&nmesh, ++numv);
			if(!v[0]->boundary()||!v[1]->boundary())
			{ 
				pt = mesh.vertexHalfedge(v[0], v[1])->he_next()->vertex()->point();
				pb = mesh.vertexHalfedge(v[1], v[0])->he_next()->vertex()->point();
				nv01->point() = (v[0]->point() + v[1]->point()) * 3 / 8 + (pt + pb) * 1 / 8;
			}
			else
			{
				nv01->point() = (v[0]->point() + v[1]->point()) / 2;
			}
			e->string() = std::to_string(numv);
		}
		else
		{
			nv01 = nmesh.idVertex(std::stoi(e->string()));
		}

		e = mesh.vertexEdge(v[1], v[2]);
		if (e->string().empty())
		{
			nv12 = d.createVertex(&nmesh, ++numv);
			if (!v[1]->boundary() || !v[2]->boundary())
			{
				pt = mesh.vertexHalfedge(v[1], v[2])->he_next()->vertex()->point();
				pb = mesh.vertexHalfedge(v[2], v[1])->he_next()->vertex()->point();
				nv12->point() = (v[1]->point() + v[2]->point()) * 3 / 8 + (pt + pb) * 1 / 8;
			}
			else
			{
				nv12->point() = (v[1]->point() + v[2]->point()) / 2;
			}
			e->string() = std::to_string(numv);
		}
		else
		{
			nv12 = nmesh.idVertex(std::stoi(e->string()));
		}

		e = mesh.vertexEdge(v[2], v[0]);
		if (e->string().empty())
		{
			nv20 = d.createVertex(&nmesh, ++numv);
			if (!v[2]->boundary() || !v[0]->boundary())
			{
				pt = mesh.vertexHalfedge(v[2], v[0])->he_next()->vertex()->point();
				pb = mesh.vertexHalfedge(v[0], v[2])->he_next()->vertex()->point();
				nv20->point() = (v[2]->point() + v[0]->point()) * 3 / 8 + (pt + pb) * 1 / 8;
			}
			else
			{
				nv20->point() = (v[2]->point() + v[0]->point()) / 2;
			}
			e->string() = std::to_string(numv);
		}
		else
		{
			nv20 = nmesh.idVertex(std::stoi(e->string()));
		}


		int nfv[3];
		nfv[0] = v[0]->id();
		nfv[1] = nv01->id();
		nfv[2] = nv20->id();
		d.createFace(&nmesh, nfv, numf++);
		nfv[0] = nv01->id();
		nfv[1] = v[1]->id();
		nfv[2] = nv12->id();
		d.createFace(&nmesh, nfv, numf++);
		nfv[0] = nv20->id();
		nfv[1] = nv12->id();
		nfv[2] = v[2]->id();
		d.createFace(&nmesh, nfv, numf++);
		nfv[0] = nv01->id();
		nfv[1] = nv12->id();
		nfv[2] = nv20->id();
		d.createFace(&nmesh, nfv, numf++);
	}

	
	// Write out the resultant obj file
	int vObjID = 1;
	std::map<int, int> vidToObjID;

	std::ofstream os(argv[2]);

	SolidVertexIterator iter(&nmesh);

	for(; !iter.end(); ++iter)
	{
		Vertex *v = *iter;
		Point p = v->point();
		os << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
		vidToObjID[v->id()] = vObjID++;
	}
	os << "# " << (unsigned int)nmesh.numVertices() << " vertices" << std::endl;

	float u = 0.0, v = 0.0;
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
	os << "# " << (unsigned int)nmesh.numVertices() << " texture coordinates" << std::endl;

	SolidFaceIterator fiter(&nmesh);
	for(; !fiter.end(); ++fiter)
	{
		Face *f = *fiter;
		FaceVertexIterator viter(f);
		os << "f " ;
		for(; !viter.end(); ++viter)
		{
			Vertex *v = *viter;
			os << vidToObjID[v->id()] << "/" << vidToObjID[v->id()] << " ";
		}
		os << std::endl;
	}
	os.close();

	std::cout << "successfully generated file " << argv[2] << std::endl;
}
