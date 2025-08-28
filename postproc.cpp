void write_file(dmatrix &u,dmatrix &v ,dmatrix &w, dmatrix &p ,int nx,int ny,int nz,double dx,double dy,double dz);





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
