
void make_ks() ;

void make_ks()
{
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	double this_km = -1.0 ;
	for (i1 = 0 ; i1 <= MB ; i1++) {
	for (i2 = 0 ; i2 <= MB ; i2++) {
	for (i3 = 0 ; i3 <= MB ; i3++) {
	for (i4 = 0 ; i4 <= MB ; i4++) {
	for (i5 = 0 ; i5 <= MB ; i5++) {
			this_km = kp*pow(KD1,i1)*pow(KD2, i2)*pow(KD3, i3)*pow(KD4, i4)*pow(KDI,i5)*exp((double)(i1+i2+i3+i4-1)*-9.0/0.6) ;
			if (this_km > min_km)
			{
				kms[i1][i2][i3][i4][i5] = this_km ;
			}
			else
			{
				kms[i1][i2][i3][i4][i5] = 0.0 ;
			}
	}
	}
	}
	}
	}
}
