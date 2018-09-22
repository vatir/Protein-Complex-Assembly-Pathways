/*
		Program for outputting the species and reactions for a "half-proteasome" of length "size"
*/

#include <bitset>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <algorithm>

using namespace std;


// Use these MACROs to turn on both inhibitors and a single type
#define INHIBITORS
#define AlphaAlpha
//#define BetaAlpha

// Set inhibitor traits here and at the top of molecule.h
#ifdef INHIBITORS
#ifdef BetaAlpha
// Beta - Alpha Class Inhibitor
bool inhibitor_block(unsigned int i, unsigned int a, unsigned int b, int n)
{
	// If i is present then a can not be
	if ((i & a) != 0) {
		return 1;
	}
	return 0;
}

bool inhibitor_allowed(unsigned int i, unsigned int a, unsigned int b, int n)
{
	if ((b == 0)&(i == 0)) { return 1; }
	// Check that i is a subset of b, since b must be present for i to be bound
	if ((b != 0) && ((i & ~b) == 0)) {
		return 1;
	}
	return 0;
}
#endif

#ifdef AlphaAlpha
// Alpha-Alpha Class Inhibitor
bool inhibitor_block(unsigned int i, unsigned int a, unsigned int b, int n)
{
	//unsigned int shifted_a;
	//shifted_a = ((1 & a) << (n - 1)) | (a >> 1);
	// If i is present then the a to the right can not be
	if ((i & a) != 0) {
		return 1;
	}
	return 0;
}

bool inhibitor_allowed(unsigned int i, unsigned int a, unsigned int b, int n)
{
	if (i == 0) { return 1; }

	unsigned int shifted_a;
	shifted_a = ((1 & a) << (n - 1)) | (a >> 1);
	unsigned int shifted_i;
	shifted_i = ((1 & i) << (n - 1)) | (i >> 1);
	//if (((i&~shifted_a) == 0)&&(i&shifted_i == 0)&&(a!=0)) {
	if (((i&~shifted_a) == 0) && (a != 0)) {
		return 1;
	}
	return 0;
}
#endif
#endif

#ifndef INHIBITORS
bool inhibitor_allowed(unsigned int i, unsigned int a, unsigned int b, int n) { return 1; }
bool inhibitor_block(unsigned int i, unsigned int a, unsigned int b, int n) { return 0; }
#endif

#include "molecule.h"
#include "reaction.h"
#include "ring_funcs.h"

int main (int argc, char ** argv)
{
	/*
	if (argc != 4)
		{
			cout << "./make_reactions.exe ring_size reactfile speciesfile" << endl ; // input the desired ring size, the output file for the reactions, and the output file for the species
			exit (1) ;
		}

	int size = atoi(argv[1]) ; // ring length, i.e. "n".	Both the a and the b rings are this length.
	ofstream out ; // ofstream for the reaction file
	out.open (argv[2]) ;
	ofstream outS ; // ofstream for the species file
	outS.open(argv[3]) ;
	*/

	// Hard Coded Inputs
	int size = 7;
	ofstream out; // ofstream for the reaction file
	out.open("reacts");
	ofstream outS; // ofstream for the species file
	outS.open("species");

	
	vector<vector<molecule> > species ; // vector containing all species
	/*	First dimension:	size class -1.	So monomers are all in species[0], dimers all in species[1], etc. */
	/*	Second dimension:	all distinct species of that size. */
	int current_name_index = 0;

	species.resize(3*size) ; // total number of size classes is equal to 2*n, i.e. max (na) + max (nb) + max (ni)
	molecule t1(size) ; // Two temporary molecules to create our monomers
	molecule t2(size) ;

	t1.make_monomer(0) ; // make the two types of monomer: a and b
	t2.make_monomer(1) ;
	
	add(t1, species, 0, current_name_index); // add these monomers to the current species list, they are monomers, so their size class is "0"
	add(t2, species, 0, current_name_index);

	#ifdef INHIBITORS
	molecule t3(size);
	t3.make_monomer(2); // make loose inhibitor
	add(t3, species, 0, current_name_index);

	//molecule t4(size);
	//t4.make_monomer(2); // make loose inhibitor
	//t4.inhibitor = 5;
	//add(t4, species, 0, current_name_index);

	//molecule t5(size);
	//t5.make_monomer(2); // make loose inhibitor
	//t5.inhibitor = 21;
	//add(t5, species, 0, current_name_index);
	//for (int inhibitor_offset = 0; inhibitor_offset < size; inhibitor_offset++) {
	//	shift_by_1(t3);
	//	add(t3, species, 0, current_name_index);
	//}
	#endif

	vector<reaction> all_reacts ; // vector containing all the reactions
	vector<reaction> p_reacts ; // temporary vector to hold all of the reactions between two molecules
	int m2_min ; // minimum for "m2" in the sum below 
	

	/* There are two size classes:	s1 and s2 */
	/* s1 will represent the monomers, i.e. s1 = 0 */
	/* s2 ranges from 0 (monomers) to 2*size-2, i.e. to molecules that are only 1 monomer away from forming the fully complete stacked ring structure */
	/* The outer loop ranges over s2 values */
	for (int s2 = 0; s2 < 3 * size; s2++)
	{
		/* for any size class, we now react it with the two monomers */
		int s1 = 0;
		for (int m1 = 0; m1 < species[s1].size(); m1++)
		{
			for (int m2 = 0; m2 < species[s2].size(); m2++)
			{
				/* m2 is a member of size class 2 */
				/* create the productive reactions for these two molecules.	Will poplulate species vector as new molecules are created */
				p_reacts = productive_reacts(species[s1][m1], species[s2][m2], species, s1, s2, current_name_index);
			}
		}
	}
	
	/* Note that all possible species can be formed by the addition of a monomer to some other species.	So the above loop creates all species that could form in the system */
	/* This loop has to be completed before the next loop, since we cannot react the a monomer, say, with species that require sequential additions of a and b monomers until
	after these species have been created */

	/*	Report on the number of species and reactions */
	cout << "Number of unique species:" << "\t" << current_name_index << endl;
	/* Print out species */
	for (int i = 0; i < species.size(); i++)
	{
		for (int j = 0; j < species[i].size(); j++)
		{
			species[i][j].print(outS);
			// Generates and prints the inhibitor target vectors.
			//species[i][j].set_inhibitor_targets();
			//species[i][j].print_inhibitor_targets(outS);
			outS << endl;
		}
	}

	/* Now determine the reactions */
	/* s1 goes over all size classes from 0 to size - 1 */
	/* note that this covers all possible unique reactions */
	for (int s1 = 0 ; s1 < 3*size ; s1++)
		{
			/* m1 is a molecule from size class s1 */
			for (int m1 = 0 ; m1 < species[s1].size() ; m1++)
	{	
		/* s2 ranges over all possible size classes that have not yet been considered, but would still be productive */
		for (int s2 = s1 ; s2+s1 < 3*size ; s2++)
			{
				/* The conditional statements below cover the difference between s1 == s2 and s1 < s2 */
				if (s2 == s1)
		{
			m2_min = m1 ;
		}
				else
		{
		m2_min = 0 ;
		}
				for(int m2 = m2_min ; m2 < species[s2].size() ; m2++)
		{
			p_reacts.resize(0) ;
			/* collect productive reactions between m1 and m2 */
			p_reacts = productive_reacts(species[s1][m1],species[s2][m2], species, s1, s2, current_name_index) ;
			/* put those reactions into all_reacts */
			for (int i = 0 ; i< p_reacts.size(); i++)
				{
					all_reacts.push_back(p_reacts[i]) ;
				}
		}
			}
	}
		}
	/* Print out reactions */
	for (int i = 0 ; i < all_reacts.size() ; i++)
		{
			all_reacts[i].print_molec(out) ;
		}

	cout << "Number of unique reactions:" << "\t" << all_reacts.size() << endl ;

}
