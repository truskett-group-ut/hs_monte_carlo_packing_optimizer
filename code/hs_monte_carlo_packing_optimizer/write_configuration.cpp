void writeclusterconfiguration(double *xc, double *yc, double *zc)
{
	int j;
	double x_com, y_com, z_com;
	//MTRand randomparticle(4567);

	//cout << endl;
	//cout<<"Writing out initial cluster configuration..." << endl;
	//cout << endl;

	ofstream output("Cluster_Configuration.xyz", ios::app);
	output << n + 1 << endl;
	output << " " << "Written by Ryan" << endl;

	x_com = 0.0;
	for (j = 0; j<n; j++) //First set of particles
	{
		output << "  " << "C       " << *(xc + j) << "       " << *(yc + j) << "       " << *(zc + j) << endl;

		x_com = x_com + *(xc + j);
		y_com = y_com + *(yc + j);
		z_com = z_com + *(zc + j);
	}

	x_com = x_com / ((double)n);
	y_com = y_com / ((double)n);
	z_com = z_com / ((double)n);

	//output << "  " << "Ar       " << 0.0 << "       " << 0.0 << "       " << 0.0 << endl;
	output << "  " << "Ar       " << x_com << "       " << y_com << "       " << z_com << endl;

	output.close();
}