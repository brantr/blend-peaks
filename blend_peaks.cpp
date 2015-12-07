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


void blend_peaks(vector<shock> *bs, vector<tracer> *bt,vector<shock> s, vector<tracer> t, double rmax)
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
  vector<tracer> tsearch;    //buffer tracers


  vector<int> interactions;
  vector<int> n_append;
  vector<int> n_merge;
  vector<int> n_expand;

  vector<tracer>::iterator it;

  //search tree
  kdtree2 *bs_tree;
  array2dfloat bs_data;
  vector<float> xc(3);  //three dimensional position

  //search results
  kdtree2_result_vector res;

  //peak trees
  kdtree2 *peak_tree;
  array2dfloat peak_tree_data;


  //flag indicating whether a search
  //tree was built
  int flag_tree_build = 0;


  printf("Blending peaks (bs %ld s %ld)..\n",bs->size(),s.size());

  //if this is our first set of nonzero peaks
  //then add them to the peak list
  if(bs->size()==0)
  {
  	printf("Empty blended peak list.\n");

    //let's make sure they're sorted by
    //density, then id

    for(ss=0;ss<s.size();ss++)
    {
      for(tt=0;tt<s[ss].l;tt++)
        tbuf.push_back(t[s[ss].o+tt]);
      
      std::sort(tbuf.begin(), tbuf.end(), tracer_density_and_id_sort);

      //add to blended tracers
      for(tt=0;tt<tbuf.size();tt++)
      {
        //set peak index
        tbuf[tt].peak_index = tbuf[0].id;

        //add to blended tracers
        bt->push_back(tbuf[tt]);      

        //printf("adding ss %ld tt %ld d %e id %ld pid %ld\n",ss,tt,tbuf[tt].d,tbuf[tt].id,tbuf[tt].peak_index);
      }

      //exit(-1);

      //reset peak index
      s[ss].id = tbuf[0].id;

      //add peak to peak list
      bs->push_back(s[ss]);

      //destroy tbuf
      vector<tracer>().swap(tbuf);
    }
/*
  	for(ss=0;ss<s.size();ss++)
  		bs->push_back(s[ss]);
  	for(tt=0;tt<t.size();tt++)
  	  bt->push_back(t[tt]);
*/
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

        //reset shock properties
        for(tt=0;tt<s[ss].l;tt++)
          tbuf.push_back(t[s[ss].o+tt]);
      
        std::sort(tbuf.begin(), tbuf.end(), tracer_density_and_id_sort);

        //add to blended tracers
        for(tt=0;tt<tbuf.size();tt++)
        {
          //set peak index
          tbuf[tt].peak_index = tbuf[0].id;
        }
        //reset shock 
        s[ss].id = tbuf[0].id;

        //append to the new list
        sappend.push_back(s[ss]);
        for(tt=0;tt<s[ss].l;tt++)
        {
          tappend.push_back(tbuf[tt]);
          //printf("ss %ld tt %ld d %e id %ld peak_index %ld\n",ss,tt,tbuf[tt].d,tbuf[tt].id,tbuf[tt].peak_index);
        }

        printf("Added ss %ld (id = %ld) to the append list.\n",ss,s[ss].id);

        //destroy tbuf
        vector<tracer>().swap(tbuf);
      }

      //ok, there is an interaction
      if(interactions.size()>0)
      {

        //if there is only one interaction, the shock
        //just needs to be expanded.  Let's get some info

        //note for *one* interaction, the peak with the lower
        //density threshold has to be the same as the previously
        //identified peak.  We don't need to be any more careful.
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
          for(tt=0;tt<(*bs)[i].l;tt++)
            tbuf.push_back((*bt)[(*bs)[i].o+tt]);

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
          {
            for(int k=0;k<3;k++)
            {
              if(tbuf[tt].x[k]<sbuf.min[k])
                sbuf.min[k] = tbuf[tt].x[k];
              if(tbuf[tt].x[k]>sbuf.max[k])
                sbuf.max[k] = tbuf[tt].x[k];
            }

            //fix peak index
            tbuf[tt].peak_index = sbuf.id;
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
//#error Add > 1 interactions
          printf("ss %ld ninteractions %ld\n",ss,interactions.size());

          for(long ti=0;ti<interactions.size();ti++)
          {
            i = interactions[ti];
            printf("ti %ld int %ld l %10ld o %10ld d %e id %10ld\n",ti,i,(*bs)[i].l,(*bs)[i].o,(*bs)[i].d,(*bs)[i].id);
          }

          //this could take a while
          //need to do search on particles
          //to find the density cross

          //if we haven't built a tree yet, let's
          //do so
          if(!flag_tree_build)
          {
            printf("Building bs tree...\n");

            //resize data for tree build
            flag_tree_build = 1;
            bs_data.resize(extents[(*bt).size()][3]);

            //load all particle data into a tree
            for(tt=0;tt<(*bt).size();tt++)
              for(int k=0;k<3;k++)
                bs_data[tt][k] = (*bt)[tt].x[k];

            //build the tree
            bs_tree = new kdtree2(bs_data, true);

            printf("done building tree...\n");
          }

          //ok, bs_tree now contains all the particles
          //in the previously identified peaks, so now
          //we can begin the peak merging

          //add this object to the merge list
          for(tt=0;tt<s[ss].l;tt++)
          {
            tbuf.push_back(t[s[ss].o+tt]);
            tbuf[tt].peak_index = -1; //
          }

          //sort by density, then id
          std::sort(tbuf.begin(), tbuf.end(), tracer_density_and_id_sort);

          //reset group id
          s[ss].id = tbuf[0].id;

          for(long ti=0;ti<interactions.size();ti++)
          {
            i = interactions[ti];

/*            
            for(tt=0;tt<(*bs)[i].l;tt++)
              tsearch.push_back((*bt)[(*bs)[i].o+tt]);
            it = std::search(tsearch.begin(), tsearch.end(), tbuf.begin(), tbuf.begin()+1, tracer_unique);

            if(it!=tsearch.end())
              printf("tbuf was found in i %ld at %ld; tbuf.id %ld ts id %ld\n",i,(it-tsearch.begin()),tbuf[0].id,tsearch[(it-tsearch.begin())].id);

            //destroy tsearch
            vector<tracer>().swap(tsearch);
*/
            //search for each peak
            tsearch.push_back((*bt)[(*bs)[i].o]);
            it = std::search(tbuf.begin(), tbuf.end(), tsearch.begin(), tsearch.begin()+1, tracer_unique);


            //Here, we need to consider separately
            //whether the interactions are merges,
            //or if the bounding boxes have just overlapped
            if(it!=tbuf.end())
            {
              //the lower threshold peak contains the higher
              //density peak.
              printf("interaction peak %6ld (id=%10ld) is in ss = %10ld (id=%10ld)\n",i,(*bs)[i].id,ss,s[ss].id);

              for(tt=0;tt<tbuf.size();tt++)
              {
                for(int k=0;k<3;k++)
                  xc[k] = tbuf[tt].x[k];
                //find any A particles within bsq
                bs_tree->n_nearest(xc,2,res);

                if(res[0].dis<1.0e-15)
                {
                  printf("ss %ld tt %ld k %d SAME tbuf.id %ld bt.id %ld\n",ss,tt,0,tbuf[tt].id,(*bt)[res[0].idx].id);
                }else{
                  printf("ss %ld tt %ld k %d DIFFERENT tbuf.id %ld bt.id %ld\n",ss,tt,0,tbuf[tt].id,(*bt)[res[0].idx].id);
                }
#error for same, check id; otherwise build second tree and then adjust peak_index to assign

                /*
                for(long k=0;k<res.size();k++)
                  printf("ss %ld tt %ld k %ld res.dis() %e\n",ss,tt,k,res[k].dis);
                */
              }
            }else{
              //here the bounding boxes have just overlapped, and
              //the lower threshold peak does not contain the peak
              //of the higher threshold peak.
              printf("interaction peak %6ld (id=%10ld) is *NOT* in ss = %10ld (id=%10ld)",i,(*bs)[i].id,ss,s[ss].id);
              exit(-1);
            }
            //destroy tsearch
            vector<tracer>().swap(tsearch);
          }

          /*
          for(tt=0;tt<s[ss].l;tt++)
          {
            printf("tracer test tt %ld d %e id %ld\n",tt,tbuf[tt].d,tbuf[tt].id);
          }
          printf("XXXX tbuf[0].id %ld s.id %ld\n",tbuf[0].id,s[ss].id);
          */

          //it = std::find(tbuf.begin(), tbuf.end(), tracer_unique);
          /*for(tt=0;tt<s[ss].l;tt++)
          {
            if(tbuf[tt].id==s[ss].id)
            {
              printf("it's here.\n");
              break;
            }
          }*/


                    
          //destroy tbuf
          vector<tracer>().swap(tbuf);
        }

      }

      //destroy the interaction list
      vector<int>().swap(interactions);
    }

    //vacate bs and bt
    (*bs).clear();
    bt->clear();

    //bs shocks with int=1 and int>1 can
    //appear in both categories

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

  //if necessary, free memory
  //associated with tree
  if(flag_tree_build)
  {
    bs_data.resize(extents[0][0]);
    free(bs_tree);
  }
}