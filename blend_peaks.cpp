#include <stdio.h>
#include "blend_peaks.hpp"

void blend_peaks(vector<shock> *bs, vector<tracer> *bt,vector<shock> s, vector<tracer> t)
{
  long ss;
  long tt;
  printf("Blending peaks..\n");

  if(bs->size()==0)
  {
  	printf("Empty blended peak list.\n");

  	for(ss=0;tt<s.size();ss++)
  		bs->push_back(s[ss]);
  	for(tt=0;tt<t.size();tt++)
  		bt->push_back(t[tt]);
  }
}