#ifndef HS_SIMULATION_H
#define HS_SIMULATION_H

class HardSphereSimulation
{
private:
	int N; //number of particles
	double eta; //volume fraction
	double L; //box edge length
	vector< array<double, 4> > particles; //coordinates of the particles

public:

	void SetDate(int year, int month, int day);

	int getYear() { return m_year; }
	int getMonth() { return m_month; }
	int getDay()  { return m_day; }
};