/* Functions used in the ring ODE-generating code
*/

void transfer(molecule&, molecule&, int); // transfers first molecule into second with offset 
void shift_by_1(molecule&); // shifts the a and b rings of molecule "m" by 1 
int isomorphic(molecule&, molecule&);  // returns "1" if two molecules are isomorphic up to offset, 0 otherwise
int isomorphic(reaction&, reaction&); // returns "1" if two reactions are isomorphic up to the order of reactants, 0 otherwise
int equiv(molecule&, molecule&);  // returns "1" if two molecules are identical, 0 otherwise

vector<reaction> productive_reacts(molecule&, molecule&, vector<vector<molecule> >&, int, int, int&);  // produces a list of productive reactions between the molecules
void react(molecule&, molecule&, molecule&);  // reacts first two molecules to form the third
//int n_gaps(molecule&); // returns the number of "gaps" in a structure.  CURRENTLY BROKEN
void add(molecule&, vector<vector<molecule> >&, int, int&); // adds the molecule to the current vector of species, and gives it a name based on the int& species index

int embeddings(molecule&, molecule&); // number of embeddings of molecule 1 into molecule 2
int embeds(molecule&, molecule&);// returns 1 if molecule 1 embeds into molecule 2, 0 otherwise
void calc_rev_coeff(reaction&); // yet another try
void subtract(molecule&, molecule&, molecule&); // subtracts molecule 1 from molecule 2, puts the result in molecule 3

// transfer takes a molecule "source" and copies it into the molecule "dest" with offset "offset"
// the offset is added to the index of every element of the vectors in "source" to obtain the corresponding element in "dest"
// so, if source.a = (1,0,0) and offset = 1, then we have dest.a = (0,1,0)
void transfer(molecule& source, molecule& dest, int offset)
{
	/* get ring length for source */
	int n = source.n;
	/* first, copy source to dest.  Takes care of resizing and other things */
	dest.copy(source);
	/* did the call request an offset greater than the ring length? */
	if (offset > n)
	{
		/* This probably should  not happen, so warn the user about it */
		cout << "WARNING:  offset greater than vector size!  Using addition on the modular integer addition group...." << endl;
		/* it is nonetheless allowed; we account for it by using addition modulo n: */
		offset %= source.n;
	}
	for (int i = 1; i <= offset; i++)
	{
		shift_by_1(dest);
	}
}

void shift_by_1(molecule& m)
{
	m.inhibitor = ((m.inhibitor & 1) << m.n - 1) | (m.inhibitor >> 1);
	m.a = ((m.a & 1) << m.n - 1) | (m.a >> 1);
	m.b = ((m.b & 1) << m.n - 1) | (m.b >> 1);
}

/* isomorphic(molecule, molecule) tells us if two molecules, "m1" and "m2", are equivalent up to rotational permutation (i.e. equivalent up to "offset") */
/* returns "1" if they are, "0" if they are not */
int isomorphic(molecule& m1, molecule& m2)
{
	/* first, we know that they CANNOT be isomorphic if they do not have the same size or number of a and b proteins */
	if (((m1.s != m2.s) || (m1.na != m2.na)) || (m1.nb != m1.nb) || (m1.ni != m1.ni))
	{
		return 0;
	}
	/* If the sizes, na and nb are identical, then they could be isomorphic, so we check if they are equivalent at 0 offset */
	else if (equiv(m1, m2))
	{
		return 1;
	}
	/* If they are not, then we check all other offsets */
	else
	{
		molecule tmp;
		tmp.copy(m2);
		for (int o = 1; o < m1.n; o++)
		{
			shift_by_1(tmp);
			if (equiv(m1, tmp))
			{
				return 1;
			}
		}
		return 0;
	}
}

/* iso_can(molecule, molecule) tells us if two molecules are isomorphic based on their canonical names */
/* WILL NOT WORK IF THE TWO MOLECULES ARE NOT ALREADY CANONICALLY NAMED */
bool iso_can(molecule& m1, molecule& m2)
{
	return (m1.name == m2.name);
}

/* isomorphic(reaction, reaction) tells us if two reactions are isomorphic up to the order of the reactants */
/* returns 1 if they are isomorphic, 0 otherwise */
/* in this implementation, we rely on the NAMES of the reactants, rather than their structures, to perform the comparison,
since this reduces the number of function calls to isomorphic(molecule, molecule) */
int isomorphic(reaction& react1, reaction& react2)
{
	if (react1.p.name != react2.p.name) { return 0; }
	if ((react1.r1.name == react2.r1.name) && (react1.r2.name == react2.r2.name))
	{
		return 1;
	}
	else if ((react1.r1.name == react2.r2.name) && (react1.r2.name == react2.r1.name))
	{
		return 1;
	}
	return 0;
}

/* equiv tells us if two molecules "m1" and "m2" are absolutely identical */
/* returns "1" if they are, "0" if they are not */
int equiv(molecule& m1, molecule&m2)
{
	return ((m1.a == m2.a) && (m1.b == m2.b) && (m1.inhibitor == m2.inhibitor));
}

/* productive_reacts is the main function of the algorithm */
/* takes in two molecules, "m1" and "m2" and determines all of the productive ways in which they can react */
/* also uses the "species" vector to determine which species is represented by the product of the reaction, or if a new species is created */
/* "s1" and "s2" are the size classes of the molcules "m1" and "m2", respectively */
/* current_name_index is the current species naming index, and is passed by reference in case we need to "add" new species */
/* returns a vector whose elements are of the class "reaction".  These are all of the productive reactions between m1 and m2 */
vector<reaction> productive_reacts(molecule& m1, molecule& m2, vector<vector<molecule> >& species, int s1, int s2, int& current_name_index)
{
	vector<reaction> res; // result: vector of reactions 
	molecule tmp; // temporary molecule used for checking all offsets of "m2"
	molecule product; // temporary molecule used for storing the product of the reaction between "m1" and "tmp"
	int can_react = 0; // variable that is used to indicate if "m1" and "tmp" can react at all
	int new_size_class = s1 + s2 + 1; // size class of the newly created molecule.  Since s1 = m1.s -1, and s2 = m2.s -1, we have s1+s2 = m1.s + m2.s -2, and we compensate by adding "1"
	int created = -1; // variable that is used to indicate if a new species is created by the reaction

	int naa_formed = 0; // number of "a-a" bonds formed by a prospective reaction
	int nbb_formed = 0; // number of "b-b" bonds formed by a prospective reaction
	int nab1_formed = 0; // number of "a-b" bonds (b to the right of a) formed by a prospective reaction
	int nab2_formed = 0; // number of "b-a" bonds (b to the left of a) formed by a prospective reaction
	int ni_formed = 0; // number of inhibitor formed by a prospective reaction
	int n_tot_formed = 0; // total number of bonds formed by a prospective reaction.  used to check if the reaction is actually productive, replaces "n_gaps" in the current implementation

	int react_already_found = 0;
	int r1_embed = -1;
	int r1_auto = -1;
	double r1_ratio = 0;
	int r2_embed = -1;
	int r2_auto = -1;
	double r2_ratio = 0;

	/* first, copy m2 into tmp */
	tmp.copy(m2);
	/* apply all possible offsets to "m2" */
	for (int offset = 0; offset < m1.n; offset++)
	{
		/* temporary reaction used to store the reaction in the current offset */
		reaction tmpr;
		/* "can_react set" to 1, will be modified to 0 if we find these molecules cannot react */
		can_react = 1;
		/* if we are not looking at the 0 offset, we shift tmp by 1 */
		if (offset != 0)
		{
			shift_by_1(tmp);
		}
		/* now we check to see if there are clashes */
		if (!(tmp.a&m1.a) && !(tmp.b&m1.b) && !(tmp.inhibitor&m1.inhibitor) && (!inhibitor_block(m1.inhibitor, tmp.a, tmp.b, tmp.n)) && (!inhibitor_block(tmp.inhibitor, m1.a, m1.b, m1.n)))
		{
			/* there are no clashes, so we can try to react them! */
			/* perform the reaction and put the result in "product" */
			react(m1, tmp, product);
			/* calculate the number of each type of bond that is formed */
			naa_formed = product.naa_bonds() - tmp.naa_bonds() - m1.naa_bonds();
			nbb_formed = product.nbb_bonds() - tmp.nbb_bonds() - m1.nbb_bonds();
			nab1_formed = product.nab1_bonds() - tmp.nab1_bonds() - m1.nab1_bonds();
			nab2_formed = product.nab2_bonds() - tmp.nab2_bonds() - m1.nab2_bonds();
			ni_formed = product.nib_bonds() - tmp.nib_bonds() - m1.nib_bonds();
			/* calculate the total number of bonds that are formed */
			n_tot_formed = naa_formed + nbb_formed + nab1_formed + nab2_formed + ni_formed;
			/* if no bonds are formed, than this reaction is not productive, since nothing actually happened */
			/* otherwise, we go ahead */
			if ((n_tot_formed != 0) && inhibitor_allowed(product.inhibitor, product.a, product.b, product.n))
			{
				/* store the two reactant molecules in the correct locations in tmpr */
				tmpr.r1 = m1;
				tmpr.r2 = m2;
				product.canonical();

				/* determine if the product actually exists already */
				/* the "created" variable pefroms two functions */
				/* if the molecule exists, then it stores the index for this product in the species[new_size_class] vector of molecules */
				/* if it does not exist, created remains "-1" which flags the need to create a new species in the species vector */
				created = -1;
				/* "j" ranges over all existing molecules in "new_size_class" */
				for (int j = 0; j < species[new_size_class].size(); j++)
				{
					/* if the product is isomorphic to some existing species, rember that species' value of j */
					if (iso_can(product, species[new_size_class][j]))
					{
						created = j;
						break;
					}
				}
				/* if we need to create a new molecule, then do so */
				if (created == -1)
				{
					/* add the product molecule to "species" in size class new_size_class, using "current_name_index" to make the new name */
					add(product, species, new_size_class, current_name_index);
					/* now, the index of the new molecule is the last in the species[new_class_size] vector, and we store that in "created" */
					created = species[new_size_class].size() - 1;
				}
				/* assign the product of the reaction to tmpr.p */
				tmpr.p = species[new_size_class][created];
				tmpr.forward_coefficient = 0;

				/* now we check to see if we have already observed this particular reaction */
				react_already_found = 0;
				for (int q = 0; q < res.size(); q++)
				{
					if (isomorphic(tmpr, res[q]))
					{
						/* if it has been, then we know there is another offset that gives this same reaction, so the forward coefficient should increase by 1 */
						react_already_found = 1;
						res[q].forward_coefficient++;
					}
				}
				if (!react_already_found)
				{
					/* If we have not already found this reaction, then this is the first offset for this case */
					tmpr.forward_coefficient = 1;
					/* put the number of bonds of each type into the correct container in tmpr */
					tmpr.bonds.push_back(naa_formed);
					tmpr.bonds.push_back(nbb_formed);
					tmpr.bonds.push_back(nab1_formed);
					tmpr.bonds.push_back(nab2_formed);
					tmpr.bonds.push_back(ni_formed);
					/* we now calculate the reverse coefficient for this reaction */
					calc_rev_coeff(tmpr);
					res.push_back(tmpr);
				}
			}
		}
	}
	return (res);
}

/* react takes two reactant molecules, "r1" and "r2", and reacts them, placing the product in "product" */
void react(molecule& r1, molecule& r2, molecule& product)
{
	product.reinit(r1.n); // initialize the product molecule with the correct sizes
						  /* if you missed a clash, shame on you! */

	if ((r1.a&r2.a) || (r1.b&r2.b) || (r1.inhibitor&r2.inhibitor))
	{
		cout << "Attempted to react two molecules with a clash:\t" << r1.name << " and " << r2.name << endl;
		exit(1);
	}
	product.inhibitor = r1.inhibitor | r2.inhibitor;
	product.a = r1.a | r2.a;
	product.b = r1.b | r2.b;
	/* calculate the product's size */
	product.ni = r1.ni + r2.ni;
	product.na = r1.na + r2.na;
	product.nb = r1.nb + r2.nb;
	product.s = product.na + product.nb + product.ni;
}

/* add takes a molecule "m" and adds it to the species vector in the size class "sc_index".  "current_name_index" is used to generate a new, unique name for the molecule */
/* the name has the convention: "Sx", where x is the number */
void add(molecule& m, vector<vector<molecule> >& species, int sc_index, int& current_name_index)
{
	/* first, increase current_name_index for this new molecule */
	current_name_index++;
	/* now put this bad boy in canonical representation */
	m.canonical();
	/* calculate the number of bonds for this structure, since it will be used frequently */
	int tmp;
	tmp = m.naa_bonds();
	tmp = m.nbb_bonds();
	tmp = m.nab1_bonds();
	tmp = m.nab2_bonds();
	tmp = m.nib_bonds();

	/* add the newly named molecule to the vector */
	species[sc_index].push_back(m);
}

/* embeddings returns the total number of ways to embed molecule "m1" into molecule "m2" */
int embeddings(molecule& m1, molecule& m2)
{
	int res = 0;
	molecule tmp;
	tmp.copy(m1);
	/* loop over all offsets of m1 */
	for (int o = 0; o < m1.n; o++)
	{
		/* if we are not at offset "0", then shift by 1 */
		if (o)
		{
			shift_by_1(tmp);
		}
		if (embeds(tmp, m2))
		{
			/* if this offset embeds in m2, then we have a valid embedding */
			res++;
		}
	}

	return res;
}

/* embeds returns 1 if molecule "m1", at its current offset, embeds into molecule "m2", and 0 otherwise */
/* a valid embedding is defined by the fact that there are no sites in "m1" that are occupied that are
not occupied in m2.  We can say that "m1 \subeq m2" if it embeds */
int embeds(molecule& m1, molecule& m2)
{
	return (((m1.a&m2.a) == m1.a) && ((m1.b&m2.b) == m1.b) && ((m1.inhibitor&m2.inhibitor) == m1.inhibitor));
}

/* calc_rev_coeff calculates the coefficient for the reverse of the given reaction. */
/* this is essentially the number of ways that the product "reaction.p" can fall
appart to give the two reactants "reaction.r1" and "reaction.r2" */
void calc_rev_coeff(reaction& reaction)
{
	molecule tmp1; // tmp1 will hold reactant "r2" at some offset "o"
	molecule tmp2; // tmp2 will hold the result of subtractiong "tmp1" from the product "p"
	molecule Orig_r1; // Copy of the original r1 (corrected for AlphaAlpha bond misalignment)
	molecule Orig_r2; // Copy of the original r2 (corrected for AlphaAlpha bond misalignment)
	int n = 0; // holds the number of valid reversals 
	int auto_r2 = 0; // automorphisms of reactant "r2"
	int equiv_r1_r2 = 0; // holds the number of embeddings of r2 into r1 
	int AlphaLocation = pow(2, (reaction.r2.n - 2)); // Numerical Location of the single Alpha inhibitor

	Orig_r1.copy(reaction.r1);
	Orig_r1.canonical();
	Orig_r2.copy(reaction.r2);
	Orig_r2.canonical();

	#ifdef AlphaAlpha
		if (((Orig_r1.inhibitor & AlphaLocation) == AlphaLocation)&&(Orig_r1.a == 0)) {
			shift_by_1(Orig_r1);
			Orig_r1.canonical();
		}
		if (((Orig_r2.inhibitor & AlphaLocation) == AlphaLocation) && (Orig_r2.a == 0)) {
			shift_by_1(Orig_r2);
			Orig_r2.canonical();
		}
	#endif

	tmp1.copy(Orig_r2);

	/* loop over possible offsets */
	for (int o = 0; o < reaction.r2.n; o++)
	{

		/* if we are not at offset "0", shift by 1 */
		if (o)
		{
			shift_by_1(tmp1);
		}
		/* see if "tmp1" will embed in "p" -- that is, see if we can "find" reactant "r2" inside "p" at offset "o" */
		/* NOTE: in our construction, "r2" is always a larger molecule than "r1".  This makes it harder for "r2" to
		embed in "p".  Since the choice of reactant here is arbitrary, trying to subtract "r2" generally reduces
		the number of computations needed */
		if (embeds(tmp1, reaction.p))
		{
			/* if we can, then we "reverse" the reaction by removing "tmp1" from "p", creating "tmp2" */
			subtract(tmp1, reaction.p, tmp2);
			tmp2.canonical();

			/* if "tmp2" is isomorphic to "r1", then we have successfully reversed the reaction */
			if (iso_can(tmp2, Orig_r1))
			{
				n++;
			}
		}
	}
	/* Now, if "r2" has automorphisms (i.e. it can embed into itself), then this will artificially inflate the
	reverse reaction coefficient.  This is because some subset of offsets will generate the same exact "tmp2"
	and thus the same exact version of the reversal.  So we normalize by the automorphisms to account for this fact */
	auto_r2 = embeddings(Orig_r2, Orig_r2);
	/* Now, if r2 is isomorphic to r1, we have to divide by an additional factor of 2 to compensate */
	if (iso_can(Orig_r1, Orig_r2))
	{
		equiv_r1_r2 = 1;
	}
	reaction.reverse_coefficient = double(n) / double(auto_r2*(1 + equiv_r1_r2));
	if (reaction.reverse_coefficient != int(reaction.reverse_coefficient)) {
		cout << "Embeddings are not integers: Something is probably wrong!!!" << "\n";
	}
}

/* subtract removes molecule "m1" from "m2" and places the resulting molecule in "res" */
void subtract(molecule& m1, molecule& m2, molecule& res)
{
	/* First make "res" a copy of "m2" */
	res.copy(m2);
	if (!embeds(m1, m2))
	{
		cout << "Attempted subtraction on a non-embedded molecule:\t" << m1.name << " - " << m2.name << endl;
		exit(1);
	}
	res.inhibitor = m2.inhibitor^m1.inhibitor;
	res.a = m2.a^m1.a;
	res.b = m2.b^m1.b;
	/* we update the total numbers of a and b proteins after the subraction */
	res.ni = m2.ni - m1.ni;
	res.na = m2.na - m1.na;
	res.nb = m2.nb - m1.nb;
	res.s = res.na + res.nb + res.ni;
}
