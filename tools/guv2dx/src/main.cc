#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <omp.h>
using namespace std;

#define MAGICNUMBER 3
#define TYPEOS 200

int main(int argc, char * argv[]) {
  string set_fname(char[]);

  int ch;
  float xs, ys, zs;
  
  while ((ch = getopt(argc, argv, "x:y:z:")) != -1) {
    switch (ch){
    case 'x':
      xs = atof(optarg);
      break;
    case 'y':
      ys = atof(optarg);
      break;
    case 'z':
      zs = atof(optarg);
      break;
    }
  }
  if (argc == 1) {
    cout << "No parameter file!" << endl ;
    return (1) ;
  }

  char* filename = argv[optind];
  std::string fname = set_fname(filename);
  ifstream fin(filename, ios::in | ios::binary);

  int nv;
  int g[3];
  double b[3];
  double s[3];
  fin.read((char *) &nv, sizeof(int));
  fin.read((char *) &g[0], sizeof(int)*3);
  fin.read((char *) &b[0], sizeof(double)*3);
  fin.read((char *) &s[0], sizeof(double)*3);

  int grid = g[0] * g[1] * g[2];
  double * guv = new double[grid * nv];
  double * tmp = new double[grid];
  fin.read((char *) &guv[0], sizeof(double) * grid * nv);
  fin.close();

  float xmin = - b[0] / 2 + xs - s[0];
  float ymin = - b[1] / 2 + ys - s[1];
  float zmin = - b[2] / 2 + zs - s[2];

  for (int iv = 0; iv < nv; ++iv) {
    char ndx[10];
    sprintf(ndx, "-%d.dx", iv);
    ofstream fout((fname + ndx).c_str(), ios::out);
    fout << "object 1 class gridpositions counts     "
	 << g[0] << "     " << g[1] << "     " << g[2] << endl;
    fout << "origin    "
	 << xmin << "   " << ymin << "   " << zmin << endl;
    fout << "delta       " << b[0] / g[0] << " 0 0" << endl;
    fout << "delta  0      " << b[1] / g[1] << " 0" << endl;
    fout << "delta  0 0      " << b[2] / g[2] << endl;
    fout << "object 2 class gridconnections counts     "
	 << g[0] << "     " << g[1] << "     " << g[2] << endl;
    fout << "object 3 class array type double rank 0 items "
	 << grid << " data follows" << endl;
    
    int ip = grid * iv;
#pragma omp parallel for
    for (int ig = 0; ig < grid; ++ig) {
      int ix = ig % g[0];
      int iy = (ig / g[0]) % g[1];
      int iz = ig / (g[0] * g[1]);
      int it = iz + iy * g[2] + ix * g[1] * g[2];
      tmp[it] = guv[ip + ig];
    }    
    for (int ig = 0; ig < grid; ++ig) {
      fout << tmp[ig] << endl;
    }
    fout << "object \"Untitled\" call field" << endl;
    fout.close();
  }
  return 0;
}
