#include <stdio.h>
#include <algorithm>
#include "blend_peaks.hpp"
#include "box_collision.hpp"

bool tracer_unique(tracer a, tracer b)
{
  return (a.id==b.id);
}
bool tracer_id_sort(tracer a, tracer b)
{
  return (a.id<b.id);
}
bool tracer_density_and_id_sort(tracer a, tracer b)
{
  if(a.d==b.d)
    return (a.id<b.id);
  return (a.d>b.d);
}
void blend_peaks(vector<shock> *bs, vector<tracer> *bt,vector<shock> s, vector<tracer> t)
{
  long ss;
  long tt;
  long ssb;
  long ttb;
  long i,j;
  long nslim;

  vector<shock>  sappend; //new shocks to add to the list
  vector<tracer> tappend;

  vector<shock>  sexpand; //shocks that need to be expanded
  vector<tracer> texpand;

  vector<shock>  smerge; //shocks that need to be merged
  vector<tracer> tmerge;

  shock  sbuf;            //buffer shocks
  vector<tracer> tbuf;    //buffer tracers

  vector<int> interactions;
  vector<int> n_append;
  vector<int> n_merge;
  vector<int> n_expand;

  vector<tracer>::iterator it;


  printf("Blending peaks (bs %ld s %ld)..\n",bs->size(),s.size());

  //if this is our first set of nonzero peaks
  //then add them to the peak list
  if(bs->size()==0)
  {
  	printf("Empty blended peak list.\n");

  	for(ss=0;ss<s.size();ss++)
  		bs->push_back(s[ss]);
  	for(tt=0;tt<t.size();tt++)
  	  bt->push_back(t[tt]);
  }else{

    //this is our main work loop

    //begin loop over peaks
    for(ss=0;ss<s.size();ss++)
      printf("s[%ld].l %ld\n",ss,s[ss].l);

//    for(ss=0;ss<s.size();ss++)
    nslim = 2;
    if(s.size()<nslim)
      nslim = s.size();
    for(ss=0;ss<nslim;ss++)
    {
      printf("Shock %ld ***************\n",ss);
      //compare bounding boxes
      for(ssb=0;ssb<bs->size();ssb++)
      {
        if(box_collision(s[ss].min,s[ss].max,(*bs)[ssb].min,(*bs)[ssb].max))
          interactions.push_back(ssb);
      }
      printf("Shock = %ld (l=%ld), number of interactions = %ld\n",ss,s[ss].l,interactions.size());

      //if there are no interactions, this is a new peak
      //this is the easiest case
      if(interactions.size() ==0)
      {
        printf("Appending shock %ld.\n",ss);
        sappend.push_back(s[ss]);
        for(tt=0;tt<s[ss].l;tt++)
          tappend.push_back(t[s[ss].o+tt]);

        printf("Added ss %ld to the append list.\n",ss);
      }

      //ok, there is an interaction
      if(interactions.size()>0)
      {

        //if there is only one interaction, the shock
        //just needs to be expanded.  Let's get some info
        if(interactions.size()==1)
        {
          i = interactions[0];
          printf("Blended shock box %e %e %e %e %e %e\n",(*bs)[i].min[0],(*bs)[i].min[1],(*bs)[i].min[2],(*bs)[i].max[0],(*bs)[i].max[1],(*bs)[i].max[2]);
          printf("New shock box     %e %e %e %e %e %e\n",s[ss].min[0],s[ss].min[1],s[ss].min[2],s[ss].max[0],s[ss].max[1],s[ss].max[2]);
          printf("Blended shock properties l %10ld o %10ld d %e id %10ld\n",(*bs)[i].l,(*bs)[i].o,(*bs)[i].d,(*bs)[i].id);
          printf("New shock properties     l %10ld o %10ld d %e id %10ld\n",s[ss].l,s[ss].o,s[ss].d,s[ss].id);

          //add this object to the merge list
          for(tt=0;tt<s[ss].l;tt++)
            tbuf.push_back(t[s[ss].o+tt]);
          for(tt=0;tt<(*bs)[ss].l;tt++)
            tbuf.push_back((*bt)[(*bs)[ss].o+tt]);

          //retain only unique tracers
          std::sort( tbuf.begin(), tbuf.end(), tracer_id_sort);
          it = std::unique(tbuf.begin(), tbuf.end(), tracer_unique);
          tbuf.resize( std::distance(tbuf.begin(), it) );

          //sort by density, then id
          std::sort(tbuf.begin(), tbuf.end(), tracer_density_and_id_sort);

          //create the shock
          sbuf.l  = tbuf.size();
          sbuf.o  = tmerge.size();
          sbuf.d  = tbuf[0].d;
          sbuf.id = tbuf[0].id;

          sbuf.min[0] = 1.e9;
          sbuf.min[1] = 1.e9;
          sbuf.min[2] = 1.e9;
          sbuf.max[0] = -1.e9;
          sbuf.max[1] = -1.e9;
          sbuf.max[2] = -1.e9;

          for(tt=0;tt<tbuf.size();tt++)
            for(int k=0;k<3;k++)
            {
              if(tbuf[tt].x[k]<sbuf.min[k])
                sbuf.min[k] = tbuf[tt].x[k];
              if(tbuf[tt].x[k]>sbuf.max[k])
                sbuf.max[k] = tbuf[tt].x[k];
            }

          //append tbuf to tmerge
          it = tmerge.end();
          tmerge.insert(it,tbuf.begin(),tbuf.end());

          smerge.push_back(sbuf);

          //destroy tbuf
          vector<tracer>().swap(tbuf);

        }


        if(interactions.size()>1)
        {
          //we have more than one peak merging
          //together
#error Add > 1 interactions
        }

      }

      //destroy the interaction list
      vector<int>().swap(interactions);
    }

    //vacate bs and bt
    (*bs).clear();
    bt->clear();

    //merge grown shocks
    for(ss=0;ss<smerge.size();ss++)
      bs->push_back(smerge[ss]);
    for(tt=0;tt<tmerge.size();tt++)
      bt->push_back(tmerge[tt]);

    //append new shocks
    for(ss=0;ss<sappend.size();ss++)
      bs->push_back(sappend[ss]);
    for(tt=0;tt<tappend.size();tt++)
      bt->push_back(tappend[tt]);

    //fix offsets
    (*bs)[0].o = 0;
    for(ss=1;ss<(*bs).size();ss++)
      (*bs)[ss].o = (*bs)[ss-1].o + (*bs)[ss].l;

  }
}