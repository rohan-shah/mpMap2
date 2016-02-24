#include "getFunnel.h"
void getFunnel(long line, Rcpp::IntegerVector& mother, Rcpp::IntegerVector& father, int* funnel, int nFounders)
{
	int currentLine = line;
	while(father(currentLine) == mother(currentLine))
	{
		currentLine = father(currentLine)-1;
	}
	int motherRow = mother(currentLine), fatherRow = father(currentLine);
	if(nFounders == 2)
	{
		funnel[0] = motherRow;
		funnel[1] = fatherRow;
		return;
	}
	int mmRow = mother(motherRow-1), mfRow = father(motherRow-1), fmRow = mother(fatherRow-1), ffRow = father(fatherRow-1);
	if(nFounders == 4)
	{
		funnel[0] = mmRow;
		funnel[1] = mfRow;
		funnel[2] = fmRow;
		funnel[3] = ffRow;
		return;
	}
	if(nFounders >= 8)
	{
		int mmmRow = mother(mmRow-1), mmfRow = father(mmRow-1), mfmRow = mother(mfRow-1), mffRow = father(mfRow-1);
		int fmmRow = mother(fmRow-1), fmfRow = father(fmRow-1), ffmRow = mother(ffRow-1), fffRow = father(ffRow-1);
		funnel[0] = mmmRow;
		funnel[1] = mmfRow;
		funnel[2] = mfmRow;
		funnel[3] = mffRow;
		funnel[4] = fmmRow;
		funnel[5] = fmfRow;
		funnel[6] = ffmRow;
		funnel[7] = fffRow;
		if(nFounders == 16)
		{
		  int mmmmRow = mother(mmmRow-1), mmmfRow=father(mmmRow-1), mmfmRow=mother(mmfRow-1), mmffRow=father(mmfRow-1);
		  int mfmmRow=mother(mfmRow-1), mfmfRow=father(mfmRow-1), mffmRow=mother(mffRow-1), mfffRow=father(mffRow-1);
		  int fmmmRow=mother(fmmRow-1), fmmfRow=father(fmmRow-1), fmfmRow=mother(fmfRow-1), fmffRow=father(fmfRow-1);
		  int ffmmRow=mother(ffmRow-1), ffmfRow=father(ffmRow-1), fffmRow=mother(fffRow-1), ffffRow=father(fffRow-1);
		  funnel[0] = mmmmRow;
		  funnel[1] = mmmfRow;
		  funnel[2] = mmfmRow;
		  funnel[3] = mmffRow;
		  funnel[4] = mfmmRow;
		  funnel[5] = mfmfRow;
		  funnel[6] = mffmRow;
		  funnel[7] = mfffRow;
		  funnel[8] = fmmmRow;
		  funnel[9] = fmmfRow;
		  funnel[10] = fmfmRow;
		  funnel[11] = fmffRow;
		  funnel[12] = ffmmRow;
		  funnel[13] = ffmfRow;
		  funnel[14] = fffmRow;
		  funnel[15] = ffffRow;
		}
	}
	else
	{
		throw std::runtime_error("Input nFounders of getFunnel must be 2, 4, 8 or 16");
	}
}
