void read_params(FILE*) ;

void read_params(FILE * in)
{
	int nin = 0 ;
	char buf[BUFSIZ] ;
	char param[20] ;
	char value[20] ;
	
	while (fgets(buf,sizeof buf,in) != NULL)
    {
		
		sscanf(buf, "%s %s", param, value) ;
		
		if(!strcmp("KD1",param))
		{
			KD1 = atof(value) ;
			nin++ ;
		}
		if(!strcmp("KD2",param))
		{
			KD2 = atof(value) ;
			nin++ ;
		}
		if(!strcmp("KD3",param))
		{
			KD3 = atof(value) ;
			nin++ ;
		}
		if(!strcmp("KD4",param))
		{
			KD4 = atof(value) ;
			nin++ ;
		}
		if(!strcmp("KDI",param))
		{
			KDI = atof(value) ;
			nin++ ;
		}
		if(!strcmp("kp",param))
		{
			kp = atof(value) ;
			nin++ ;
		}
		if(!strcmp("kpi",param))
		{
			kpi = atof(value) ;
			nin++ ;
		}
		if(!strcmp("min_km",param))
		{
			min_km = atof(value) ;
			nin++ ;
		}
		if(!strcmp("tstart",param))
		{
			tstart = atof(value) ;
			nin++ ;
		}
		if(!strcmp("tfac",param))
		{
			tfac = atof(value) ;
			nin++ ;
		}
		if(!strcmp("A0",param))
		{
			A0 = atof(value) ;
			nin++ ;
		}
		if(!strcmp("B0",param))
		{
			B0 = atof(value) ;
			nin++ ;
		}
		if(!strcmp("I0",param))
		{
			I0 = atof(value) ;
			nin++ ;
		}
		if(!strcmp(param,"nout"))
		{
			nout = atoi(value) ;
			nin++ ;
		}
		if(!strcmp(param,"n" ))
		{
			n = atoi(value) ;
			nin++ ;
		}
		if(!strcmp(param,"abstol"))
		{
			ATOL = atof(value) ;
			nin++ ;
		}
		if(!strcmp(param,"reltol"))
		{
			RTOL = atof(value) ;
			nin++ ;
		}
		if(!strcmp(param,"mstep"))
		{
			maxstep = atoi(value) ;
			nin++ ;
		}
    }
	if (nin != INPARAMS)
    {
		printf("%s\t%d\n","Incorrect number of parrameters input", nin) ;
		exit(1) ;
    }
}
