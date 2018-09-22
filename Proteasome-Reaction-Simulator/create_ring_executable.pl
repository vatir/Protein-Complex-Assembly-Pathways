#!/usr/bin/perl

if (@ARGV != 5)
  {
    die "./parse_mathematican n species_file reaction_file name_of_outstruct unique_name\n" ;
  }
$n = $ARGV[0] ;
open (INS, "$ARGV[1]") ;
open (INR, "$ARGV[2]") ;
$nameout = $ARGV[3] ;
$code_root = "2_stacked_rings_n_".$n."_".$ARGV[4] ;
$code_filename = $code_root .".c" ;
open (OUT, ">$code_filename") ;

@species = () ;
@ode_terms = () ;
@init_conds = () ;
%spec_map = () ;

$max_bond = 0 ;
$outstruct_ind = -1 ;
$i = 0 ;
while (<INS>)
  {
    chomp ;
    @l = split(/;;/) ;
    $name = $l[0] ;
    if ($name == $nameout)
      {
	$outstruct_ind = $i ;
      }
    $spec_map{$name} = $i ;
    $species[$i] = $name ;
    $i++ ;
  }
$tot_reacts = $i ;
print "$tot_reacts\n" ;

while (<INR>)
  {
    chomp ;
    @l = split(/\t/); 
    @r1code = split(/;;/,$l[0]) ;
    $r1 = $r1code[0] ;
    $r1_index = $spec_map{$r1} ;
    @r2code = split(/;;/,$l[1]) ;
    $r2 = $r2code[0] ;
    $r2_index = $spec_map{$r2} ;
    $forward_coeff = $l[2] ;
    @pcode = split(/;;/,$l[3]) ;
    $p = $pcode[0] ;
    $p_index = $spec_map{$p} ;
    @bonds = split(/,/,$l[4]) ;
    for($i = 0 ; $i < @bonds ; $i++)
      {
	if ($bonds[$i] > $max_bond)
	  {
	    $max_bond = $bonds[$i] ;
	  }
      }
    $rev_coeff = $l[5] ;
    if ($r1 eq $r2)
      {
	$forward_div = "0.5" ;
      }
    else
      {
	$forward_div = "1.0" ;
      }
    $forward_div = "1.0" ; # "Fix" factor of 2, which violates conservation of mass
    if ($bonds[4] eq 0)
	{
    	$f_term = $forward_coeff.".0 * kp*".$forward_div."*ydat[".$r1_index."]*ydat[".$r2_index."] " ;
	}
    else
	{
    	$f_term = $forward_coeff.".0 * kpi*".$forward_div."*ydat[".$r1_index."]*ydat[".$r2_index."] " ;
	}
    $r_term = $rev_coeff.".0 *kms[".$bonds[0]."][".$bonds[1]."][".$bonds[2]."][".$bonds[3]."][".$bonds[4]."]*ydat[".$p_index."] " ;
    $r1_term = "- ".$f_term."+ ".$r_term ;
    $r2_term = $r1_term ;
    $p_term = "+ ".$f_term."- ".$r_term ;
    $ode_terms[$r1_index] .= $r1_term ;
    $ode_terms[$r2_index] .= $r2_term ;
    $ode_terms[$p_index] .= $p_term ;
  }

for ($i = 0 ; $i < @species ; $i++)
  {
    $ode_terms[$i] .= " ;" ;
  }

print OUT "/* Automatically Generated file for integrating ring ODEs using CVODE */\n\n";
print OUT "#include <stdio.h>\n#include <stdlib.h>\n#include <math.h>\n#include \"llnltyps.h\"\n#include \"cvode.h\"\n#include \"cvdense.h\"\n" ;
print OUT "#include \"cvband.h\"\n#include \"cvdiag.h\"\n#include \"nvector.h\"\n#include \"llnlmath.h\"\n\n\n" ;
print OUT "/* Problem Constants */\n#define ITOL  SS\n#define ERRFP stdout\n#define OPTIN TRUE\n#define INPARAMS 18\n\n\n" ;

print OUT "#define FULLIND $outstruct_ind\n" ;
print OUT "#define MB $max_bond\n" ;
print OUT "#define NEQ $tot_reacts\n\n" ;

print OUT "double KD1 = -1.0 ;\ndouble KD2 = -1.0 ;\ndouble KD3 = -1.0 ;\ndouble KD4 = -1.0 ;\ndouble KDI = -1.0;\ndouble kp = -1.0 ;\ndouble kpi = -1.0;\ndouble min_km = -1.0 ;\ndouble tstart = -1.0 ;\n" ;
print OUT "double tfac = -1.0 ;\ndouble A0 = -1.0 ;\ndouble B0 = -1.0 ;\ndouble I0 = -1.0 ;\nint nout = -1 ;\nint n = -1  ;\ndouble ATOL = -1.0 ;\ndouble RTOL = -1.0 ;\n" ;
print OUT "int maxstep = -1.0 ;\ndouble kms[MB+1][MB+1][MB+1][MB+1][MB+1] ;\n#include \"read_params.h\"\n#include \"make_ks.h\"\n\n\n" ;

print OUT "/* Private Helper Functions */\nstatic int  Problem1(FILE*);\n\n/* Functions Called by the CVODE Solver */\n" ;
print OUT "static void f(integer N, real t, N_Vector y, N_Vector ydot, void *f_data);\n" ;

print OUT "/* Implementation */\n\nint main(int argc, char** argv)\n{\n  int nerr ;\n  double k = 0.0 ;\n  if (argc != 3)\n    {\n";
print OUT "      printf(\"%s\\n\", \"./run.exe in_params outfile\") ;\n      exit(1) ;\n    }\n  FILE* in = fopen(argv[1],\"r\") ;\n" ;
print OUT "  FILE* f_o = fopen(argv[2],\"w\") ;\n" ;
print OUT "\n  read_params(in) ;\n  make_ks() ;\n  \n  nerr = Problem1(f_o) ;\n\n  return(0);\n}\n\n\n" ;

print OUT "static int Problem1(FILE* out)\n{\n  real ropt[OPT_SIZE], reltol=RTOL, abstol=ATOL, t, tout, ero, er;\n  long int iopt[OPT_SIZE];\n" ;
print OUT "  int lmm, miter, flag, iter, iout, nerr=0;\n  N_Vector y;\n  void *cvode_mem;\n  iopt[MAXORD] = 0 ;\n  iopt[MXSTEP] = maxstep ;\n  iopt[MXHNIL] = 0 ;\n" ;
print OUT "  ropt[H0] = 0.0 ;\n  ropt[HMAX] = 0.0 ;\n  ropt[HMIN] = 0.0 ;\n\n  double outfrac = 0.0 ;\n\n  int converge = 1 ;\n\n" ;
print OUT "  y = N_VNew(NEQ, NULL);\n  N_VIth(y,0) = A0 ;\n N_VIth(y,2) = I0 ;\n  N_VIth(y,1) = B0 ;\n\n  cvode_mem = CVodeMalloc(NEQ, f, 0, y, BDF, NEWTON, ITOL,\n" ;
print OUT "			  &reltol, &abstol, NULL, ERRFP, OPTIN, iopt, ropt, NULL);\n  if (cvode_mem == NULL) { printf(\"CVodeMalloc failed.\"); return(1); }\n";
print OUT "  CVDense(cvode_mem,NULL,NULL) ;\n\n  for(iout=1, tout=tstart; iout <= nout; iout++, tout *= tfac) \n" ;
print OUT "    {\n      flag = CVode(cvode_mem, tout, y, &t, NORMAL);\n      if (flag != SUCCESS) {\n	nerr++;\n" ;
print OUT "	printf(\"\\n\\n CVode returned error flag = %d at time %g\\n\\n\", flag,t);\n	converge = 0 ;\n	break;	\n      }\n";
print OUT "      outfrac = (double)n*2.0*N_VIth(y, FULLIND)/A0 ;\n      fprintf(out, \"%g\t%g\\n\", tout,outfrac) ;\n    }\n  CVodeFree(cvode_mem);\n  N_VFree(y);\n" ;
print OUT "  return(nerr);\n}\n\n\n" ;

print OUT "static void f(integer N, real t, N_Vector y, N_Vector ydot, void *f_data)\n{\n\n" ;
print OUT "  real *ydat, *dydat ;\n\n" ;
print OUT "  ydat = N_VDATA(y) ;\n" ;
print OUT "  dydat = N_VDATA(ydot) ;\n\n" ;

for ($i = 0 ; $i < @species ; $i++)
  {
    print OUT "  dydat[$i] = $ode_terms[$i]\n" ;
  }
print OUT "}\n" ;


open (OUTM, ">Makefile") ;

print OUTM "# Makefile for a 2 stacked ring ODE integrator using CVODE\n" ;
print OUTM "# Generated automatically\n\n\n" ;

print OUTM "COMPILER = gcc\n\nINC = ./CVODE/include\n\nLIB = ./CVODE/lib\n\nOPTS = -I\$(INC)\n";
print OUTM "\nHDRS = \$(INC)/llnltyps.h \$(INC)/cvode.h \$(INC)/cvdense.h \$(INC)/cvband.h \\\n       \$(INC)/cvdiag.h \$(INC)/nvector.h \$(INC)/llnlmath.h\n" ;
print OUTM "\n".$code_root.": ".$code_root.".o\n	\$(COMPILER) -L\$(LIB) -o ".$code_root." ".$code_root.".o -lcvode -lm\n" ;
print OUTM "\n".$code_root.".o: ".$code_root.".c \$(HDRS)\n	\$(COMPILER) \$(OPTS) -c ".$code_root.".c\n" ;

#system("make") ;
