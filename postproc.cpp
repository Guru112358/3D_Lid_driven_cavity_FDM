void write_file(dmatrix &u,dmatrix &v ,dmatrix &w, dmatrix &p ,int nx,int ny,int nz,double dx,double dy,double dz);
void write_file_vtk_ascii(dmatrix &u, dmatrix &v, dmatrix &w, dmatrix &p,int nx, int ny, int nz, double dx, double dy, double dz);




void write_file(dmatrix &u,dmatrix &v ,dmatrix &w, dmatrix &p ,int nx,int ny,int nz,double dx,double dy,double dz)
{
 	 std::string str="final";
	std::string ext=".csv";

	//str=str.append(std::to_string(number));
	str=str.append(ext);

	std::fstream data;

	data.open(str,std::ios::out);
	data<<" x , y , z, u , v , w, p"<<std::endl;
			


for (int i=0; i<nx; i++)
{
	for (int j=0; j<ny; j++)
	
	{
	for (int k=0; k<nz; k++)
		{
        	data<<i*dx<<" , "<<j*dy<<" , "<<k*dz<<" , " <<u(i,j,k)<<" , "<<v(i,j,k)<<" , "<<w(i,j,k)<<" , "<<p(i,j,k)<<std::endl;
		}				
	}
	
}

    data.close();

}


void write_file_vtk_ascii(dmatrix &u, dmatrix &v, dmatrix &w, dmatrix &p,int nx, int ny, int nz, double dx, double dy, double dz)
{
    std::ofstream out("final.vtk");
    if (!out.is_open())
    {
        std::cerr << "Error: could not open output file.\n";
        return;
    }

    // --- VTK header ---
    out << "# vtk DataFile Version 3.0\n";
    out << "3D CFD results\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING " << dx << " " << dy << " " << dz << "\n";
    out << "POINT_DATA " << nx * ny * nz << "\n";

    // --- Write pressure (scalar) ---
    out << "SCALARS p double\n";
    out << "LOOKUP_TABLE default\n";
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                out << p(i, j, k) << "\n";

    // --- Write velocity (vector) ---
    out << "VECTORS velocity double\n";
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                out << u(i, j, k) << " " << v(i, j, k) << " " << w(i, j, k) << "\n";

    out.close();
}
