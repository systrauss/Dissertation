/* Quick macro to convert root histogram to ascii*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

/* ======================================================================= */
int file_error(char *error_type, char *filename)
{
  /* write error message */
  /* cannot perform operation error_type on file filename */

  if (strlen(error_type) + strlen(filename) > 58) {
    warn1("ERROR - cannot %s file\n%s\n", error_type, filename);
  } else {
    warn1("ERROR - cannot %s file %s\n", error_type, filename);
  }
  return 0;
} /* file_error */

/* ======================================================================= */
int put_file_rec(FILE *fd, void *data, int numbytes)
{
  /* write one fortran-unformatted style binary record into data */
  /* returns 1 for error */

#ifdef VMS  /* vms */
  int   j1;
  short rh[2];
  char  *buf;

  buf = data;
  j1 = numbytes;
  if (numbytes <= 2042) {
    rh[0] = numbytes + 2; rh[1] = 3;
  } else {
    rh[0] = 2044; rh[1] = 1;
    while (j1 > 2042) {
      if (fwrite(rh, 2, 2, fd) != 2 ||
	  fwrite(buf, 2042, 1, fd) != 1) return 1;
       rh[1] = 0; j1 -= 2042; buf += 2042;
    }
    rh[0] = j1 + 2; rh[1] = 2;
  }
  if (fwrite(rh, 2, 2, fd) != 2 ||
      fwrite(buf, j1, 1, fd) != 1) return 1;
  /* if numbytes is odd, write an extra (padding) byte */
  if (2*(numbytes>>1) != numbytes) {
    j1 = 0;
    fwrite(&j1, 1, 1, fd);
  }
    
#else /* unix */

  if (fwrite(&numbytes, 4, 1, fd) != 1 ||
      fwrite(data, numbytes, 1, fd) != 1 ||
      fwrite(&numbytes, 4, 1, fd) != 1) return 1;
#endif
  return 0;
} /*put_file_rec */

/* ======================================================================= */

int wspec(char *filnam, float *spec, int idim)
{
  /* write spectra in gf3 format
     filnam = name of file to be created and written
     spec = spectrum of length idim */

  char buf[32];
  int  j, c1 = 1, rl = 0;
  char namesp[8];
  FILE *file;

  file = fopen(filnam,"w");
  //if (!(file = open_new_file(filnam, 0))) return 1;
  //strncpy(namesp, filnam, 8);
  //if (j < 8) memset(&namesp[j], ' ', 8-j);

  /* WRITE(1) NAMESP,IDIM,1,1,1 */
  /* WRITE(1) SPEC */
#define W(a,b) { memcpy(buf + rl, a, b); rl += b; }
  W(namesp,8); W(&idim,4); W(&c1,4); W(&c1,4); W(&c1,4);
#undef W
cout << "trying to put it in the file..." << endl;
  if (put_file_rec(file, buf, rl) ||
      put_file_rec(file, spec, 4*idim)) {
    file_error("write to", filnam);
    fclose(file);
    return 1;
  }
  cout << "idk man" << endl;
  fclose(file);
  return 0;
} 

/* ======================================================================= */

/* Macro to convert root histogram to spe for radware*/

void r2s(TH1F* inName, char* fileOut)
{
	float spec[16384];
  	int   idim1;
 	int   i, numch;

	TH1F* specToConv = inName;
	int n = specToConv->GetNbinsX();

/**	if (file_out == NULL)
		printf("Sorry, but the ascii output file did not open.");
**/
	i = 0;
	while (i<=(n+1) && i < 16384)
	{
		//file_out << "\t" << specToConv->GetBinCenter(i) <<"," << "\t" << specToConv->GetBinContent(i) << endl;
		spec[i++] = specToConv->GetBinContent(i);
	}
	printf("%i lines read...\n", i);

	numch = idim1 = i;
	printf(" %i channels..\n", numch);

	wspec(Form("%s.spe",fileOut), spec, numch);
    /* tell user that the file has been converted */
    printf(" histogram ==> %s.spe, %i chs.\n", fileOut, numch);
}

/* ======================================================================= */

void r2a(TH1F* inName, char* fileOut)
{
	ofstream file_out(fileOut);
	TH1F* specToConv = inName;
	int n = specToConv->GetNbinsX();

/**	if (file_out == NULL)
		printf("Sorry, but the ascii output file did not open.");
**/
	for (int i=0; i<=(n+1); i++)
	{
		file_out << "\t" << specToConv->GetBinCenter(i) <<"," << "\t" << specToConv->GetBinContent(i) << endl;
	}

	file_out.close();
}

void fileR2S(char* fileIn, char* fileOut, int nStart, int nEnd)
{
  TFile* file = new TFile(fileIn);
  cout << fileIn << " open" << endl;
  for(int i=nStart;i<nEnd;i++)
  {
    //r2s((TH1F*)file->Get(Form("clover_raw_%d",i)),Form("%s%dl%d",fileOut,(int)i/4,i%4));
    //r2s((TH1F*)file->Get(Form("Clover_%d",i)),Form("%s%dl%d",fileOut,(int)i/4,i%4)); //Special case right now
  }
  file->Close();
}

void fileR2S(char* fileIn, char* fileOut)
{
  TFile* file = new TFile(fileIn);
  cout << fileIn << " open" << endl;
  for(int i=0;i<2;i++)
  {
    r2s((TH1F*)file->Get(Form("Clover_%d",i)),Form("%s_C%d",fileOut,i));
  }
  for(int i=1;i<7;i++)
  {
    r2s((TH1F*)file->Get(Form("SiLi%d",i)),Form("%s_S%d",fileOut,i));
  }
  file->Close();
}