/* Molecule class for deterministic half proteasomes
*/
#include <string>
#include <bitset>

using namespace std;

#ifdef INHIBITORS

#ifdef BetaAlpha
// Inhibitor Bond Count (Beta - Alpha)
int inhibitor_bond_count(unsigned int inhibitor, unsigned int a, unsigned int b, int n)
{
	int bonds = 0;
	for (int i = 0; i < n; i++)
	{
		if (((1 << i)&inhibitor) && ((1 << i)&b))
		{
			bonds++;
		}
	}
	return bonds;
}
#endif

#ifdef AlphaAlpha
// Inhibitor Bond Count (Alpha - Alpha)
int inhibitor_bond_count(unsigned int inhibitor, unsigned int a, unsigned int b, int n)
{
	unsigned int shifted;
	unsigned int inhibitor_count;
	
	shifted = ((1 & a) << (n - 1)) | (a >> 1);
	inhibitor_count = ((((a&(~(a << 1))) >> 1) | ((a << n - 1)&~a))&inhibitor);
	//inhibitor_count = (shifted&inhibitor);

	int bonds = 0;
	for (int i = 0; i < n; i++)
	{
		if (inhibitor_count&(1<<i))
		{
			bonds++;
		}
	}
	return bonds;
}
#endif
#endif

#ifndef INHIBITORS
int inhibitor_bond_count(unsigned int inhibitor, unsigned int a, unsigned int b, int n) { return 0; }
#endif

class molecule
{
public:
	/* Rings in the "canonical" offset */
	unsigned int a = 0; // stores a copy of the "a" ring but in a canonical offset.  used for printing and isomorphism
	unsigned int b = 0; // stores a copy of the "b" ring but in a canonical offset.  used for printing and isomorphism
	unsigned int inhibitor = 0; // stores a copy of the "i" ring but in a canonical offset.  used for printing and isomorphism

						/* Inhibitor Vectors : Bit representations of possible inhibitor binding locations in terms of the subunit that is displaced */
	unsigned int alpha_alpha = 0;	// between alphas
	unsigned int beta_beta = 0;		// between betas
	unsigned int alpha1_beta1 = 0;	// On the alpha between the alpha and the lower right beta
	unsigned int beta1_alpha1 = 0;	// On the beta between the alpha and the lower right beta
	unsigned int alpha1_beta7 = 0;	// On the alpha between the alpha and the lower left beta
	unsigned int beta7_alpha1 = 0;	// On the beta between the alpha and the lower left beta
	void set_inhibitor_targets(); // Creates the values of the above, based on the current values of a and b

								  /* Properties */
	int n; // size of the ring
	int ni = 0; // number of inhibitors currently in the molecule
	int na = 0; // number of alpha proteins currently in the molecule
	int nb = 0; // number of beta proteins currently in the molecule
	int s = 0; // total number of proteins--size class
	unsigned int name = 0; // name of the molecule

						   /* Member Functions */
	molecule(); // default constructor
	molecule(int); // constructor with the size of the rings given

	void reinit(int); // resets a molecule with the size of the rings given

	int naa_bonds(); // function that returns the number of "a-a" bonds in the molecule
	int nbb_bonds(); // function that returns the number of "b-b" bonds in the molecule
	int nab1_bonds(); // function that returns the number of "a-b" (b to the right of a) bonds in the molecule
	int nab2_bonds(); // function that returns the number of "b-a" (b to the left of a) bonds in the molecule
	int nib_bonds(); // function that returns the number of inhibitor bonds in the molecule

	void copy(molecule&); // turns this molecule into an exact copy of another
	void make_monomer(int); // makes this into a monomer.  If int == 0, then we do an "a" monomer.  If int == 1, then "b". 

	void canonical(); // obtains the offset that gives the canonical representation

	string format_vector(unsigned int);
	string comma_interleave(string in_string);

	void print(ofstream&);  // function that prints the canonical representation of the molecule in a specified format to the file indicated in ofstream
	void print_inhibitor_targets(ofstream&);  // function that prints the posible inhibitor targets in a specified format to the file indicated in ofstream

private:
	/* note, these numbers are initialized to -1 so we can calculate them only when needed */
	int naa = -1; // stores the number of "a-a" bonds in the molecule
	int nbb = -1; // stores the number of "b-b" bonds in the molecule
	int nab1 = -1; // stores the number of "a-b" (b to the right of a) bonds in the molecule
	int nab2 = -1; // stores the number of "b-a" (b to the left of a) bonds in the molecule
	int nib = -1; // stores the number of inhibitors/"inhibitors bound" in the molecule
};

/* Constructors championship */
/* Default Constructor */
molecule::molecule() {}

/* construct the molecule with a specified ring length "size" */
/* equivalent to the constructor above, but with n specified */
molecule::molecule(int size) { n = size; }

/* Initialize the molecule with a specified ring length "size" */
void molecule::reinit(int size)
{
	inhibitor = 0;
	a = 0; //delete any contents in the a and b rings
	b = 0;
	n = size;
	ni = 0;
	na = 0;
	nb = 0;
	s = na + nb + ni;
	name = 0;
	naa = -1;
	nbb = -1;
	nab1 = -1;
	nab2 = -1;
	nib = -1;
}

/* Create an exact copy of the source molecule */
void molecule::copy(molecule& source)
{
	inhibitor = source.inhibitor;
	a = source.a;
	b = source.b;
	n = source.n;
	ni = source.ni;
	na = source.na;
	nb = source.nb;
	s = na + nb + ni;
	name = source.name;
	/* THE VALUES OF naa ETC ARE NOT CHANGED HERE */
}

// Format a string for output to cout
string molecule::format_vector(unsigned int holder)
{
	string out;
	string bs_string;

	bitset<64> bs(holder);
	bs_string = bs.to_string();
	out = bs_string.substr(bs_string.length() - n, bs_string.length());
	return out;
}

/* naa_bonds returns the number of "a-a" bonds */
int molecule::naa_bonds()
{
	/* first checks if this has been done before.  if so, returns the private variable naa */
	if (naa != -1)
	{
		return naa;
	}
	else
	{
		int bonds = 0;
		/* create 2 copies of a */
		unsigned int ac1 = a;
		unsigned int ac2 = a;
		unsigned int a_shift = 0;
		/* a_shift is a copy of a, shifted right by 1 */
		a_shift = ((1 & ac1) << (n - 1)) | (ac2 >> 1);
		/* bond_vector now contains a "1" whenever i & i-1 are occupied in a */
		int bond_vector = a & a_shift;
		for (int i = 0; i < n; i++)
		{
			if ((1 << i)&bond_vector)
			{
				bonds++;
			}
		}
		/* store the value, and return it */
		naa = bonds;
		return naa;
	}
}

/* nbb_bonds returns the number of "b-b" bonds */
int molecule::nbb_bonds()
{
	/* first checks if this has been done before.  if so, returns the private variable nbb */
	if (nbb != -1)
	{
		return nbb;
	}
	else
	{
		int bonds = 0;
		/* create 2 copies of b */
		unsigned int bc1 = b;
		unsigned int bc2 = b;
		unsigned int b_shift = 0;
		/* b_shift is a copy of b, shifted right by 1 */
		b_shift = ((1 & bc1) << (n - 1)) | (bc2 >> 1);
		/* bond_vector now contains a "1" whenever i & i-1 are occupied in v */
		int bond_vector = b & b_shift;
		for (int i = 0; i < n; i++)
		{
			if ((1 << i)&bond_vector)
			{
				bonds++;
			}
		}
		/* store the value, and return it */
		nbb = bonds;
		return nbb;
	}
}

/* nab1_bonds returns the number of "a-b" (b on the right side of a) bonds */
int molecule::nab1_bonds()
{
	/* first checks if this has been done before.  if so, returns the private variable nab1 */
	if (nab1 != -1)
	{
		return nab1;
	}
	else
	{
		int bonds = 0;
		/* an "a-b" bond ONLY occurs when there is an "a" and a "b" protein at the same site "i" */
		/* so the bitwise & returns a vector that is 1 if both are occupied at "i", 0 otherwise */
		int bond_vector = a&b;
		for (int i = 0; i < n; i++)
		{
			if ((1 << i)&bond_vector)
			{
				bonds++;
			}
		}
		/* store the value, and return it */
		nab1 = bonds;
		return nab1;
	}
}

/* nab2_bonds returns the number of "b-a" (b on the left side of a) bonds */
int molecule::nab2_bonds()
{
	/* first checks if this has been done before.  if so, returns the private variable nab2 */
	if (nab2 != -1)
	{
		return nab2;
	}
	else
	{
		int bonds = 0;
		/* there is a "b-a" bond if there is a "b" at "i" and an "a" at "i+1" */
		/* to check this, we shift b right by 1, then perform a bitwise & */
		/* create 2 copies of b */
		unsigned int bc1 = b;
		unsigned int bc2 = b;
		unsigned int b_shift = 0;
		/* b_shift is a copy of b, shifted right by 1 */
		b_shift = ((1 & bc1) << (n - 1)) | (bc2 >> 1);
		int bond_vector = a&b_shift;
		for (int i = 0; i < n; i++)
		{
			if ((1 << i)&bond_vector)
			{
				bonds++;
			}
		}
		nab2 = bonds;
		return nab2;
	}
}

int molecule::nib_bonds()
{
	// Check if cached
	if (nib != -1)
	{
		return nib;
	}
	nib = inhibitor_bond_count(inhibitor, a, b, n);
	return nib;
}

/* Convert the molecule into a monomer */
/* "Flag" determines if this is an "a" or a "b" monomer */
void molecule::make_monomer(int flag)
{
	if (flag == 0)  // i.e. if flag == 0, then this will be an "a" monomer
	{
		a = 1;
		na = 1;
	}
	else if (flag == 1) 
	{
		b = 1;
		nb = 1;
	}
	else if (flag == 2)
	{
		inhibitor = (1 << (n - 2)); // Move Inhibitor over to allow for Alpha-Alpha class bonds
		ni = 1;
	}
	s = na + nb + ni;
}

/* Converts the current offset of the molecule into a canonical offset */
/* The canonical representation is the offset that gives the maximum integer value for "a" */
/* If all offsets of "a" are identical, then we maximize b*/
/* If all offsets of "a" and "b" are identical, then the molecule is already in the canonical representation */
/* The results are stored in a_can and b_can */
void molecule::canonical()
{
	#ifdef AlphaAlpha // For the Alpha-Alpha class inhibitors canonical form will never bind to anything, so don't move the single inhibitor
	//if ((inhibitor == (1 << n - 2))&&(a==0)&&(b==0)) {
	if ((inhibitor == (1 << (n - 2)))&&(a == 0) && (b == 0)) {
			name = inhibitor << 2 * n;
		return; 
	}
	#endif
	unsigned int max = inhibitor << 2*n | a << n | b;
	unsigned int new_name = 0;
	unsigned int max_i = inhibitor;
	unsigned int max_a = a;
	unsigned int max_b = b;
	for (int i = 1; i < n; i++)
	{
		inhibitor = ((inhibitor & 1) << (n - 1)) | (inhibitor >> 1);
		a = ((a & 1) << (n - 1)) | (a >> 1);
		b = ((b & 1) << (n - 1)) | (b >> 1);
		new_name = inhibitor << 2 * n | a << n | b;
		if (new_name > max)
		{
			max = new_name;
			max_i = inhibitor;
			max_a = a;
			max_b = b;
		}
	}
	inhibitor = max_i;
	a = max_a;
	b = max_b;
	name = max;
}

// Calculate Inhibitor Targets for the Species (Not currently in use)
void molecule::set_inhibitor_targets()
{
	// Inhibitor is bound to the first part of the name and the location in the vector is the displaced protein in the vector of the second part of the name
	alpha_alpha = (((a&(~(a << 1))) >> 1) | ((a << n - 1)&~a));

	beta_beta = (((b&(~(b << 1))) >> 1) | ((b << n - 1)&~b));

	// 1 and 1 means on the same index
	alpha1_beta1 = a&~b;

	// 1 and 7 means the bond that shifts to beta on the left so in the case of a 7 ring alpha1 would be binding to beta 7
	alpha1_beta7 = (a&(~(b >> 1))) << 1;
	alpha1_beta7 = (alpha1_beta7 | (alpha1_beta7 >> n))&(~b);
	//alpha1_beta7 = (alpha1_beta7 << 1) | (alpha1_beta7 >> (n - 1));

	beta1_alpha1 = b&~a;

	beta7_alpha1 = (((b&(~(a << 1))) >> 1) | ((b << n - 1)&~a));

	// Cleanup just in case they get used for something else later
	unsigned int mask = (1 << n) - 1;
	alpha_alpha = alpha_alpha&mask;
	beta_beta = beta_beta&mask;
	alpha1_beta1 = alpha1_beta1&mask;
	alpha1_beta7 = alpha1_beta7&mask;
	beta1_alpha1 = beta1_alpha1&mask;
	beta7_alpha1 = beta7_alpha1&mask;
}

// Output Inhibitor Targets for the Species (Not currently in use)
void molecule::print_inhibitor_targets(ofstream& out)
{
	out << ";";
	out << format_vector(alpha_alpha);
	out << ";";
	out << format_vector(beta_beta);
	out << ";";
	out << format_vector(alpha1_beta1);
	out << ";";
	out << format_vector(alpha1_beta7);
	out << ";";
	out << format_vector(beta1_alpha1);
	out << ";";
	out << format_vector(beta7_alpha1);
	out << ";";
}

/* Interleave commas into the input string, without a trailing comma*/
string molecule::comma_interleave(string in_string)
{
	string out_string;
	for (char i_char : in_string.substr(0, in_string.length() - 1))
	{
		out_string += i_char;
		out_string += ",";
	}
	out_string += in_string.substr(in_string.length() - 1, in_string.length());
	return out_string;
}

/* Print the molecule to an ofstream with the format:
name;;a[n-1],a[n-2],...,a[0];b[n-1],b[n-2],...,b[0] */
void molecule::print(ofstream& out)
{
	string i_str = format_vector(inhibitor);
	string a_str = format_vector(a);
	string b_str = format_vector(b);

	out << name << ";;";
	out << comma_interleave(i_str);
	out << ";";
	out << comma_interleave(a_str);
	out << ";";
	out << comma_interleave(b_str);
}
