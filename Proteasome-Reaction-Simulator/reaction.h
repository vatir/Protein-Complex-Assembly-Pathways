/* Molecule class for deterministic half proteasomes
*/

using namespace std ;

class reaction
{
	public:
	molecule r1 ; // reactant molecule 1
	molecule r2 ; // reactant molecule 2
	molecule p ; // product molecule
	vector<int> bonds ; // vector of the number of each type of bond formed  FORMAT:  aa, bb, ab1, ab2
	
	double reverse_coefficient ; // number of ways of performing the reverse reaction
	double forward_coefficient ; // number of ways of performing the forward reaction

	void print_names(ofstream&) ; // print out a simple species-name version of the reaction
	void print_molec(ofstream&) ; // print out the reaction using species codes
	private:
};

/* Print the reaction just with the names */
void reaction::print_names(ofstream& out)
{
	out << r1.name << "\t " << r2.name << "\t" << forward_coefficient << "\t" << p.name << "\t(" ;
	for (int i = 0 ; i < bonds.size()-1 ; i++)
	{
	out << bonds[i] << "," ;
	}
	out << bonds[bonds.size()-1] <<")\t" << reverse_coefficient << endl ;
}

/* Print the reaction using the format for molecule.print() */
/* overall format is:
	reactant1 \t reactant2 \t forward_coefficient \t product \t naa,nbb,nab1,nab2,nib \t reverse_coefficient */
/* forward_coefficient should multiply "k_plus*r1*r2" in the ODEs to take care of the number of ways r1 and r2 can react */
/* naa is the number of "a-a" bonds formed by the reaction */
/* nbb is the number of "b-b" bonds formed by the reaction */
/* nab1 is the number of "a-b" (b is on the right side of a) bonds formed by the reaction */
/* nab2 is the number of "b-a" (b is on the left side of a) bonds formed by the reaction */
/* reverse_coefficient should multiply "k_minus(naa,nbb,nab1,nab2)*p" in the ODEs to represent the number of ways to find 
	r1 and r2 in p */
void reaction::print_molec(ofstream& out)
{
	/* Print out the code for the first reactant */
	r1.print(out) ;
	out << "\t" ;
	/* Print out the code for the second reactant */
	r2.print(out) ;
	/* Print out the forward_coefficient */
	out << "\t" << forward_coefficient << "\t" ;
	p.print(out) ;
	out << "\t" ;
	
	/* Print out the bonds formed vector */
	for (int i = 0 ; i < bonds.size()-1 ; i++)
	{
	out << bonds[i] << "," ;
	}
	out << bonds[bonds.size()-1] << "\t" ;
	/* Print out the reverse_coefficient */
	out  << reverse_coefficient << endl ;
}  


