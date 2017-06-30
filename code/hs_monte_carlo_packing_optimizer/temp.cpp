//checks for overlaps based on cell list
/*bool CheckParticleOverlap(State &state, Particle &particle){
	//if call list is active use it
	if (active)
	{
		//cell that the particle resides in
		int cell_x_part, cell_y_part, cell_z_part;
		tie(cell_x_part, cell_y_part, cell_z_part) = FindCell(particle);

		//neighbor cells
		int cell_x, cell_y, cell_z;
		int cell_x_wrpd, cell_y_wrpd, cell_z_wrpd;

		//loop over neighbor cells
		for (int cell_x = cell_x_part - 1; cell_x <= cell_x_part + 1; cell_x++){
			for (int cell_y = cell_y_part - 1; cell_y <= cell_y_part + 1; cell_y++){
				for (int cell_z = cell_z_part - 1; cell_z <= cell_z_part + 1; cell_z++){
					cell_x_wrpd = cell_x + (int)(cell_x == -1)*N_cells - (int)(cell_x == N_cells)*N_cells;
					cell_y_wrpd = cell_y + (int)(cell_y == -1)*N_cells - (int)(cell_y == N_cells)*N_cells;
					cell_z_wrpd = cell_z + (int)(cell_z == -1)*N_cells - (int)(cell_z == N_cells)*N_cells;
					//loop over particles 
					for (auto it = cells[cell_x_wrpd][cell_y_wrpd][cell_z_wrpd].begin(); it != cells[cell_x_wrpd][cell_y_wrpd][cell_z_wrpd].end(); it++){

					}
				}
			}
		}
	}
	//brute force is needed as cell list is no good
	else{

	}



}*/