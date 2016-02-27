#include <stdio.h>
#include <algorithm>
#include <iostream>
#include "blend_peaks.hpp"
#include "box_collision.hpp"

struct edgeg
{
  long tt;
  long idxA;
  float disA;
  long idxg;
  float disg;
  long pid;
};
struct tracer_key
{
  long  idx;
  float peak_density;
  long  peak_index;
  long  id;
  float d;
};
bool tracer_key_sort(tracer_key a, tracer_key b)
{
  if(a.peak_density==b.peak_density)
  {
    if(a.peak_index==b.peak_index)
    {
      if(a.d==b.d)
      {
        return (a.id<b.id);

      }else{
        return (a.d > b.d);
      }
    }else{
      return (a.peak_index < b.peak_index);
    }
  }else{
    return (a.peak_density>b.peak_density);
  }
}
bool edgeg_pid_sort(edgeg a, edgeg b)
{
  return (a.pid<b.pid);
}
bool edgeg_disA_sort(edgeg a, edgeg b)
{
  return (a.disA<b.disA);
}
bool edgeg_disg_sort(edgeg a, edgeg b)
{
  return (a.disg<b.disg);
}
bool edgeg_pid_assigned(edgeg a)
{
  return (a.pid!=-1);
}
bool tracer_unique(tracer a, tracer b)
{
  return (a.id==b.id);
}
bool tracer_id_sort(tracer a, tracer b)
{
  return (a.id<b.id);
}
bool shock_id_sort(shock a, shock b)
{
  return (a.id<b.id);
}
bool shock_density_sort(shock a, shock b)
{
  return (a.d>b.d);
}
bool tracer_pid_and_density_and_id_sort(tracer a, tracer b)
{
  if(a.peak_index==b.peak_index)
  {
    if(a.d==b.d)
      return (a.id<b.id);
    return (a.d>b.d);
  }
  return (a.peak_index<b.peak_index);
}
bool tracer_density_and_id_sort(tracer a, tracer b)
{
  if(a.d==b.d)
    return (a.id<b.id);
  return (a.d>b.d);
}
bool tracer_position(tracer a, tracer b)
{
  if(a.x[0]==b.x[0])
  {
    if(a.x[1]==b.x[1])
    {
      return(a.x[2]<b.x[2]);
    }else{
      return(a.x[1]<b.x[1]);
    }
  }
  return (a.x[0]<b.x[0]);
}

void set_peak_box(shock *s, vector<tracer> t)
{
  for(int k=0;k<3;k++)
  {
    s->min[k] =  1.0e9;
    s->max[k] = -1.0e9;
  }

  for(long tt=0;tt<t.size();tt++)
    for(int k=0;k<3;k++)
    {
      if(t[tt].x[k]<s->min[k])
        s->min[k] = t[tt].x[k];
      if(t[tt].x[k]>s->max[k])
        s->max[k] = t[tt].x[k];
    }
}

void keep_duplicates(vector<tracer> tunion, vector<tracer> *toverlap)
{
  vector<tracer>::iterator ia;

  ia = std::adjacent_find(tunion.begin(), tunion.end(), tracer_unique);
  if(ia!=tunion.end())
  {
    toverlap->push_back(*ia);
    while(ia!=tunion.end())
    {
      ia = std::adjacent_find(++ia, tunion.end(), tracer_unique);
      if(ia!=tunion.end())
        toverlap->push_back(*ia);
    }
  }
}

void keep_unique(vector<tracer> tcbuf, vector<tracer> *tcorr)
{
  long tt=0;
  if(tcbuf.size()==1)
    tcorr->push_back(tcbuf[0]);
  while(tt<tcbuf.size()-1)
  {
    if(tcbuf[tt].id!=tcbuf[tt+1].id)
    {
      tcorr->push_back(tcbuf[tt]);
    }else{
      tt++;
    }
    tt++;
  }
  if(tcbuf.size()>1)
    if(tcbuf[tcbuf.size()-2].id!=tcbuf[tcbuf.size()-1].id)
      tcorr->push_back(tcbuf[tcbuf.size()-1]);
}

void blend_peaks(vector<shock> *bs, vector<tracer> *bt,vector<shock> s, vector<tracer> t, double rmax)
{
  long ss;
  long tt;
  long ssb;
  long ttb;
  long i,j;
  long nslim;

  long id_check = 56992049; //debugging

  vector<shock>  sappend; //new shocks to add to the list
  vector<tracer> tappend;

  vector<shock>  sexpand; //shocks that need to be expanded
  vector<tracer> texpand;

  vector<shock>  smerge; //shocks that need to be merged
  vector<tracer> tmerge;

  vector<shock>  sstore; //final revision before storing in bs
  vector<tracer> tstore; //final revision before storing in bt

  shock  sbuf;            //buffer shocks
  vector<tracer> tbuf;    //buffer tracers
  vector<tracer> tsearch; //buffer tracers


  vector<tracer> tunion;    //union of high and low rho threshold peaks
  vector<tracer> toverlap;  //overlap between high and low rho threshold peaks
  vector<tracer> tcorr;       //corrected high-density threshold peak
  vector<tracer> tcorr_lowd;  //corrected low density threshold peak
  vector<tracer> tcbuf;

  vector<shock>  scomp;       //compare blended shocks before and after


  vector<int> interactions;
  //vector<int> n_append;       //THESE ARE NEW LDT SHOCKS
  //vector<int> n_merge;
  //vector<int> n_expand;       //THESE ARE HDT SHOCKS WITH ADDED PARTICLES

  vector<tracer>::iterator it;
  vector<tracer>::iterator ia;

  long pcount = 0;
  long tbb = 0;
  long tba = 0;

  //search tree
  kdtree2 *bs_tree;
  array2dfloat bs_data;
  vector<float> xc(3);  //three dimensional position

  //search results
  kdtree2_result_vector res;

  //peak trees
  kdtree2 *peak_tree;
  array2dfloat peak_tree_data;

  //gap trees
  kdtree2 *gap_tree;
  array2dfloat gap_tree_data;

  vector<tracer> tgap;
  int flag_gap_tree = 0;

  //search results
  kdtree2_result_vector gap_res;

  //flag indicating whether a search
  //tree was built
  int flag_tree_build = 0;

  int flag_multi = 0; //a shock had multiple interactions

  int flag_box_fail = 0;  //a shock had multi interactions because boxes, not tracers, overlapped

  vector<long> lcheck;
  vector<long>::iterator il;

  printf("Blending peaks (bs %ld s %ld t %ld)..\n",bs->size(),s.size(),t.size());

  //verify that tracers in t are unique

  for(tt=0;tt<t.size();tt++)
  {
  	lcheck.push_back(t[tt].id);
    if(t[tt].id==id_check)
      printf("tt %ld id %ld\n",tt,t[tt].id);
  }
  std::sort(lcheck.begin(), lcheck.end());
  il = std::unique(lcheck.begin(), lcheck.end());
  lcheck.resize( std::distance(lcheck.begin(), il) );
  if(lcheck.size()!=t.size())
  {
    printf("*********** lcheck.size() %ld t.size() %ld\n",lcheck.size(),t.size());
    exit(-1);

  }

  //cin.get();

  //#error The particles are actually lost. sum s.l == t.size(), all t are unique.

  //perform a very simple check on the particle ids

  //if this is our first set of nonzero peaks
  //then add them to the blended peak list
  if(bs->size()==0)
  {
  	printf("Empty blended peak list.\n");

    //let's make sure they're sorted by
    //density, then id

    //loop over shocks 
    for(ss=0;ss<s.size();ss++)
    {
      //some simple checking
      pcount += s[ss].l;


      //add all the tracers from
      //this shock to a buffer
      for(tt=0;tt<s[ss].l;tt++)
      {
        tbuf.push_back(t[s[ss].o+tt]);

        //if(t[s[ss].o+tt].id==5708317)
        if((t[s[ss].o+tt].id==id_check)||(t[s[ss].o+tt].id==22663579)||(t[s[ss].o+tt].id==43349505))
        {
          printf("PRESENT A pid %ld ss %ld s.id %ld\n",t[s[ss].o+tt].id,ss,s[ss].id);
          //exit(-1);
        }
      }
      
      //sort all the tracers by density, and then by id
      std::sort(tbuf.begin(), tbuf.end(), tracer_density_and_id_sort);

      //create the new, blended shock
      sbuf.l  = tbuf.size();
      sbuf.o  = tappend.size();
      sbuf.d  = tbuf[0].d;
      sbuf.id = tbuf[0].id;
      
      //set the peak box
      set_peak_box(&sbuf, tbuf);

      //adjust the tracers' peak indices
      for(tt=0;tt<tbuf.size();tt++)
        tbuf[tt].peak_index = sbuf.id;

      //append tbuf to tappend
      it = tappend.end();
      tappend.insert(it,tbuf.begin(),tbuf.end());

      //append sbuf to smerge
      sappend.push_back(sbuf);

      //destroy tbuf
      vector<tracer>().swap(tbuf);

    }//end loop over shocks
  }else{//end bs->size()==0

    //this is our main work loop
    //and bs->size()>0

    //save the blended shock list for comparison
    for(ss=0;ss<bs->size();ss++)
      scomp.push_back((*bs)[ss]);

    //begin loop over peaks found at this density threshold
    for(ss=0;ss<s.size();ss++)
      printf("s[%4ld].l %6ld id %10ld\n",ss,s[ss].l,s[ss].id);


    //Begin loop over all shocks
    //found in the shock list from the
    //current density threshold

    //we have to decide what to do with
    //them based on whether they contain
    //peaks already identified at higher
    //density thresholds or not.

    nslim = s.size();
    for(ss=0;ss<nslim;ss++)
    {

      
      //reset the flag of multiple interactions
      flag_multi = 0;

      //compare bounding boxes
      for(ssb=0;ssb<bs->size();ssb++)
      {
        if(box_collision(s[ss].min,s[ss].max,(*bs)[ssb].min,(*bs)[ssb].max))
          interactions.push_back(ssb);
      }

      printf("******************Shock %ld %10ld***************\n",ss,s[ss].id);
      printf("Shock = %6ld (l=%6ld; id=%10ld; d=%e), number of interactions = %ld\n",ss,s[ss].l,s[ss].id,s[ss].d,interactions.size());
      

      //if there are no interactions, this is a new peak
      //this is the easiest case, add shock to sappend
      //and add tracers to tappend

      //For this case, there should not be any 
      //particles in *bt that need to be added to
      //the list
      if(interactions.size() ==0)
      {
      	//some simple checking
      	pcount += s[ss].l;

        //printf("Appending shock %ld.\n",ss);

        //Put tracers into a buffer
        for(tt=0;tt<s[ss].l;tt++)
        {
          tbuf.push_back(t[s[ss].o+tt]);
          /*if(t[s[ss].o+tt].id==5708317)
          {
            printf("PRESENT B ss %ld\n",ss);
            exit(-1);
          }*/
          if(t[s[ss].o+tt].id==id_check)
          {
  //          printf("PRESENT B ss %ld\n",ss);
            printf("PRESENT B pid %ld ss %ld s.id %ld\n",t[s[ss].o+tt].id,ss,s[ss].id);

          }

        }
      
        //sort tracers by density, then id
        std::sort(tbuf.begin(), tbuf.end(), tracer_density_and_id_sort);

        //set peak index
        for(tt=0;tt<tbuf.size();tt++)
          tbuf[tt].peak_index = tbuf[0].id;

        //reset shock peak index
        s[ss].id = tbuf[0].id;

        //set the peak box
        set_peak_box(&s[ss], tbuf);

        //append shock and tracers to the append lists
        sappend.push_back(s[ss]);
        for(tt=0;tt<s[ss].l;tt++)
          tappend.push_back(tbuf[tt]);

        //printf("Added shock properties ss %ld l %10ld o %10ld d %e id %10ld\n",ss,s[ss].l,s[ss].o,s[ss].d,s[ss].id);
        printf("************\n");
        printf("*** Blended (A; int==0) shock properties l %10ld o %10ld d %e id %10ld\n",s[ss].l,s[ss].o,s[ss].d,s[ss].id);
        printf("************\n");


        //destroy tracer buffer
        vector<tracer>().swap(tbuf);

      }//end interactions==0

      //ok, there is at least one interaction
      if(interactions.size()>0)
      {

        //if there is only one interaction and the peak is the same, the shock
        //just needs to be expanded.  Let's get some info

        //note for *one* interaction, the peak with the lower
        //density threshold is usually to be the same as the previously
        //identified peak.  Then, texpand and sexpand contain the lists of the
        //simply grown shocks

        if((interactions.size()==1)&&(s[ss].id==(*bs)[interactions[0]].id) )
        //if((interactions.size()==1)&&((s[ss].id==(*bs)[interactions[0]].id)||(s[ss].d==(*bs)[interactions[0]].d)) )
        {

          //we have found one of the HDT peaks as a LDT peak that
          //has no other interactions.  We can then just duplicate
          //the HDT peak and expand it

          //To do so, we have to add all the LDT particles from *t
          //to the *bt particles in the HDT peak, and then 
          //compare to find duplicates -- only keeping unique
          //particles


          //some simple checking
      	  pcount += s[ss].l;

          i = interactions[0];
          //printf("EXPANDING:\n");
          //printf("HDT shock box %e %e %e %e %e %e\n",(*bs)[i].min[0],(*bs)[i].min[1],(*bs)[i].min[2],(*bs)[i].max[0],(*bs)[i].max[1],(*bs)[i].max[2]);
          //printf("LDT shock box %e %e %e %e %e %e\n",s[ss].min[0],s[ss].min[1],s[ss].min[2],s[ss].max[0],s[ss].max[1],s[ss].max[2]);
          //printf("HDT shock properties l %10ld o %10ld d %e id %10ld\n",(*bs)[i].l,(*bs)[i].o,(*bs)[i].d,(*bs)[i].id);
          //printf("LDT shock properties l %10ld o %10ld d %e id %10ld\n",s[ss].l,s[ss].o,s[ss].d,s[ss].id);

          //add the high and low density threshold
          //shocks' tracers to a buffer
          for(tt=0;tt<s[ss].l;tt++)
          {
            tbuf.push_back(t[s[ss].o+tt]);
          }
          for(tt=0;tt<(*bs)[i].l;tt++)
            tbuf.push_back((*bt)[(*bs)[i].o+tt]);

          //retain only unique tracers to blend the peaks
          std::sort( tbuf.begin(), tbuf.end(), tracer_id_sort);
          it = std::unique(tbuf.begin(), tbuf.end(), tracer_unique);
          tbuf.resize( std::distance(tbuf.begin(), it) );

          //sort by density, then id
          std::sort(tbuf.begin(), tbuf.end(), tracer_density_and_id_sort);

          //create the new, blended shock
          sbuf.l  = tbuf.size();
          sbuf.o  = texpand.size();
          sbuf.d  = tbuf[0].d;
          sbuf.id = tbuf[0].id;


          //set the peak box
          set_peak_box(&sbuf, tbuf);

          //adjust the tracers' peak indices
          for(tt=0;tt<tbuf.size();tt++)
            tbuf[tt].peak_index = sbuf.id;

          //append tbuf to texpand
          it = texpand.end();
          texpand.insert(it,tbuf.begin(),tbuf.end());

          //append sbuf to sexpand
          sexpand.push_back(sbuf);

          //destroy tbuf
          vector<tracer>().swap(tbuf);
          printf("************\n");
          printf("*** Blended (B, int==1) shock properties l %10ld o %10ld d %e id %10ld\n",sbuf.l,sbuf.o,sbuf.d,sbuf.id);
          printf("************\n");

        //}else if((interactions.size()==1)&&(s[ss].d<(*bs)[interactions[0]].d)){
        }else if(interactions.size()==1){


#ERROR IF THE PROXIMATE INTERACTION IS a HDT WITH A LOWER DENSITY PEAK than the LDT peak, IT WILL BE SUBSUMED!

          //proximate shocks with different peaks
          //so we need to add ss as if interactions==0

          //THESE ARE REALLY THE SAME PEAK, BLENDED TOGETHER BY LOW D LIKELY

          //Nominally, this involves adding the LDT 
          //particles from the input *t.

          i = interactions[0];
          printf("PROXIMATE\n");
          printf("HDT shock box %e %e %e %e %e %e\n",(*bs)[i].min[0],(*bs)[i].min[1],(*bs)[i].min[2],(*bs)[i].max[0],(*bs)[i].max[1],(*bs)[i].max[2]);
          printf("LDT shock box %e %e %e %e %e %e\n",s[ss].min[0],s[ss].min[1],s[ss].min[2],s[ss].max[0],s[ss].max[1],s[ss].max[2]);
          printf("HDT shock properties l %10ld o %10ld d %e id %10ld\n",(*bs)[i].l,(*bs)[i].o,(*bs)[i].d,(*bs)[i].id);
          printf("LDT shock properties l %10ld o %10ld d %e id %10ld\n",s[ss].l,s[ss].o,s[ss].d,s[ss].id);

          //some simple checking
          pcount += s[ss].l;

          //Put tracers into a buffer
          for(tt=0;tt<s[ss].l;tt++)
          {
            tbuf.push_back(t[s[ss].o+tt]);

            if(t[s[ss].o+tt].id==(*bs)[i].id)
            {
              printf("HDT PEAK INDEX PRESENT IN PROXIMATE C pid %ld ss %ld s.id %ld\n",t[s[ss].o+tt].id,ss,s[ss].id);

            }
          }
      
          //sort tracers by density, then id
          std::sort(tbuf.begin(), tbuf.end(), tracer_density_and_id_sort);

          //set peak index
          for(tt=0;tt<tbuf.size();tt++)
            tbuf[tt].peak_index = tbuf[0].id;

          //reset shock peak index
          s[ss].id = tbuf[0].id;

          //set the peak box
          set_peak_box(&s[ss], tbuf);

          //append shock and tracers to the append lists
          sappend.push_back(s[ss]);
          for(tt=0;tt<s[ss].l;tt++)
            tappend.push_back(tbuf[tt]);

          //printf("Added shock properties ss %ld l %10ld o %10ld d %e id %10ld\n",ss,s[ss].l,s[ss].o,s[ss].d,s[ss].id);
          printf("************\n");
          printf("*** Blended (C, int=1) shock properties l %10ld o %10ld d %e id %10ld\n",s[ss].l,s[ss].o,s[ss].d,s[ss].id);
          printf("************\n");

          //destroy tracer buffer
          vector<tracer>().swap(tbuf);

        }

        //if a shock has more than one interaction
        //enter this conditional. If this is turned off,
        //shocks with multiple interactions will not 
        //remain
        if(interactions.size()>1)
        {
          //some simple checking
      	  pcount += s[ss].l;

          //we have more than one peak merging
          //together
          printf("ss %ld ninteractions %ld\n",ss,interactions.size());

          for(long ti=0;ti<interactions.size();ti++)
          {
            i = interactions[ti];
            printf("ti %6ld int %6ld l %10ld o %10ld d %e id %10ld\n",ti,i,(*bs)[i].l,(*bs)[i].o,(*bs)[i].d,(*bs)[i].id);
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

          //insert the tracers in the low density
          //threshold peak into a buffer called tbuf
          //and set the peak indices to -1
          for(tt=0;tt<s[ss].l;tt++)
          {
            tbuf.push_back(t[s[ss].o+tt]);
            tbuf[tt].peak_index = -1;
//            if(t[s[ss].o+tt].id==5708317)
            if((t[s[ss].o+tt].id==id_check)||(t[s[ss].o+tt].id==22663579)||(t[s[ss].o+tt].id==43349505))
            {
//              printf("PRESENT D ss %ld\n",ss);
              printf("PRESENT D pid %ld ss %ld s.id %ld\n",t[s[ss].o+tt].id,ss,s[ss].id);

              //exit(-1);
            }
          }

          //sort by density, then id
          std::sort(tbuf.begin(), tbuf.end(), tracer_density_and_id_sort);

          //create the new, blended shock
          sbuf.l  = tbuf.size();
          sbuf.o  = 0;
          sbuf.d  = tbuf[0].d;
          sbuf.id = tbuf[0].id;

          //begin a loop over all the interactions
          for(long ti=0;ti<interactions.size();ti++)
          {
            //set i to the index of the hdt peak
            //that has a box collision with the
            //ldt peak
            i = interactions[ti];

            //search for each high density threshold peak
            //in the tracers of the low density threshold peak
            tsearch.push_back((*bt)[(*bs)[i].o]);
            it = std::search(tbuf.begin(), tbuf.end(), tsearch.begin(), tsearch.begin()+1, tracer_unique);


            //WE have several possible conditions:
            //
            //We have probably merged together a
            //variety of peaks.  We want to keep
            //the dense peaks distinct, but appropriately
            //assign particles from the low density 
            //medium around the HDT peaks to one of the
            //HDT peaks

            //One possibility is to remove all the
            //HDT peak particles from this low density
            //peak, and then appropriately assign the
            //remaining particles to the HDT peaks


            //Here, we need to consider separately
            //whether the interactions are merges,
            //or if the bounding boxes have just overlapped

//CASE 1 -- we find the low density peak in a high density peak
            if((it!=tbuf.end()) && !flag_multi)
            {
              //remember that we've done this multiple interaction
              flag_multi = 1;

              //in this branch,
              //the lower threshold peak contains the higher
              //density peak.


              printf("**INT #%ld: LDT DOES    CONTAIN HDT int peak %6ld (id=%10ld) is       in ss = %10ld (id=%10ld)\n",ti,i,(*bs)[i].id,ss,s[ss].id);

              //why not just build a separate tree buffer?
              //loop over tbuf, search the tree for nearest
              //particle, which is either in a peak (then assign it)
              //or is not, and we place the particle in the 
              //array tgap
              for(tt=0;tt<tbuf.size();tt++)
              {
               if(tbuf[tt].id==id_check)
                {

                  printf("PRESENT E pid %ld ss %ld s.id %ld\n",tbuf[tt].id,ss,s[ss].id);

                  //exit(-1);
                }
                //check the location of the particle
                for(int k=0;k<3;k++)
                  xc[k] = tbuf[tt].x[k];

                //find any A particles within bsq
                bs_tree->n_nearest(xc,2,res);

                //if we've found the same particle, go ahead
                //and assign it to the HDT peak
                if((*bt)[res[0].idx].id==tbuf[tt].id)
                {


                  if(tbuf[tt].id==id_check)
                  {
                    printf("PRESENT E in BT pid %ld ss %ld s.id %ld\n",tbuf[tt].id,ss,s[ss].id);

                    //exit(-1);
                  }
                  //assign to HDT
                  tbuf[tt].peak_index = (*bt)[res[0].idx].peak_index;

                }else{
                  if(tbuf[tt].id==id_check)
                  {
                    printf("PRESENT E in TGAP PB pid %ld ss %ld s.id %ld\n",tbuf[tt].id,ss,s[ss].id);

                    //exit(-1);
                  }
                  //add this tracer to the gap region
                  tgap.push_back(tbuf[tt]);
                }
              }

              printf("tgap.size() %ld\n",tgap.size());

              //if there are particles in the gap (likely), 
              //build a tree for the particles in the gap
              if(tgap.size()>0)
              {
                edgeg ein;
                vector<edgeg> egap;
                vector<edgeg>::iterator ei;

                //build tree for tracers in gap between peaks
                //(or around peaks at least)
                gap_tree_data.resize(extents[tgap.size()][3]);
                for(tt=0;tt<tgap.size();tt++)
                {

                  if(tgap[tt].id==id_check)
                  {
                    printf("PRESENT E in gap pid %ld ss %ld s.id %ld\n",tgap[tt].id,ss,s[ss].id);

                    //exit(-1);
                  }
                  tgap[tt].peak_index = -1; //reset 
                  for(int k=0;k<3;k++)
                    gap_tree_data[tt][k] = tgap[tt].x[k];
                }
                gap_tree = new kdtree2(gap_tree_data, true);
                flag_gap_tree = 1;

                //loop over particles in the gap
                for(tt=0;tt<tgap.size();tt++)
                {
                  //initialize query vectory
                  for(int k=0;k<3;k++)
                    xc[k] = gap_tree_data[tt][k];

                  //search on HDT
                  bs_tree->n_nearest(xc,1,res);

                  //search on gap particles
                  gap_tree->n_nearest(xc,2,gap_res);

                  //if the HDT tracers are closer than
                  //any gap tracers, go ahead and assign
                  //this tracer to the HDT shock
                  if(res[0].dis<gap_res[1].dis)
                    tgap[tt].peak_index = (*bt)[res[0].idx].peak_index;


                  //if we still don't know the peak
                  //add this to the list of tracers
                  //we need to sort out
                  if(tgap[tt].peak_index==-1)
                  { 
                    ein.tt  = tt;
                    ein.idxA = res[0].idx;
                    ein.disA = res[0].dis;
                    ein.idxg = gap_res[1].idx;
                    ein.disg = gap_res[1].dis;
                    ein.pid  = -1;
                    egap.push_back(ein);
                  }
                }//tgap.size()

                printf("egap.size() %ld\n",egap.size());

                int neiter=0;

                //iteratively "bind" tracers
                //in gap to a shock
                /*slow but works
                while(egap.size()>0)
                {
                  //printf("****\n");

                  neiter++;
                  if(!(neiter%100))
                    printf("iteration %d egap %ld\n",neiter,egap.size());

                  //OK, let's continue
                  //sort the edges by distance to peaks
                  std::sort(egap.begin(),egap.end(),edgeg_disA_sort);

                  //set the peak_index of the tracer that is closest to a peak
                  tgap[egap[0].tt].peak_index = (*bt)[egap[0].idxA].peak_index;
                  egap[0].pid                 = tgap[egap[0].tt].peak_index;

                  //find the tracers that have this tracer as its
                  //nearest neighbor, and assign those to this shock
                  //as well.
                  for(tt=0;tt<egap.size();tt++)
                  {
                    if(tgap[egap[tt].idxg].peak_index!=-1)
                    { 
                      tgap[egap[tt].tt].peak_index = tgap[egap[tt].idxg].peak_index;
                      egap[tt].pid = tgap[egap[tt].idxg].peak_index;
                    }
                    //printf("tt %ld egap.tt %ld tracer id %ld d %e pid %10ld bres id %ld bres dis %e bres d %e bres pid %ld gres id %ld gres dis %e gres d %e gres pid %ld\n",tt,egap[tt].tt,tgap[egap[tt].tt].id,tgap[egap[tt].tt].d,tgap[egap[tt].tt].peak_index,(*bt)[egap[tt].idxA].id,egap[tt].disA,(*bt)[egap[tt].disA].d,(*bt)[egap[tt].disA].peak_index,tgap[egap[tt].idxg].id,egap[tt].disg,tgap[egap[tt].idxg].d,tgap[egap[tt].idxg].peak_index);
                  }

                  //let's resort and remove unassigned tracers
                  std::sort(egap.begin(),egap.end(),edgeg_pid_sort);

                  //remove tracers assigned to a peak
                  for(tt=0;tt<egap.size();tt++)
                    if(egap[tt].pid!=-1)
                     break;

                  //remove assigned tracers
                  egap.erase(egap.begin()+tt,egap.end());
                }//while(egap.size()>0)
                */

                //sort the edges by distance to peaks
                std::sort(egap.begin(),egap.end(),edgeg_disA_sort);

                while(egap.size()>0)
                {
                  //printf("****\n");

                  neiter++;
                  //if(!(neiter%100))
                    //printf("iteration %d egap %ld\n",neiter,egap.size());

                  //OK, let's continue


                  //set the peak_index of the tracer that is closest to a peak
                  tgap[egap[0].tt].peak_index = (*bt)[egap[0].idxA].peak_index;
                  egap[0].pid                 = tgap[egap[0].tt].peak_index;

                  //find the tracers that have this tracer as its
                  //nearest neighbor, and assign those to this shock
                  //as well.
                  for(tt=0;tt<egap.size();tt++)
                  {
                    if(tgap[egap[tt].idxg].peak_index!=-1)
                    { 
                      tgap[egap[tt].tt].peak_index = tgap[egap[tt].idxg].peak_index;
                      egap[tt].pid = tgap[egap[tt].idxg].peak_index;
                    }
                    //printf("tt %ld egap.tt %ld tracer id %ld d %e pid %10ld bres id %ld bres dis %e bres d %e bres pid %ld gres id %ld gres dis %e gres d %e gres pid %ld\n",tt,egap[tt].tt,tgap[egap[tt].tt].id,tgap[egap[tt].tt].d,tgap[egap[tt].tt].peak_index,(*bt)[egap[tt].idxA].id,egap[tt].disA,(*bt)[egap[tt].disA].d,(*bt)[egap[tt].disA].peak_index,tgap[egap[tt].idxg].id,egap[tt].disg,tgap[egap[tt].idxg].d,tgap[egap[tt].idxg].peak_index);
                  }

                  ei = std::remove_if(egap.begin(),egap.end(),edgeg_pid_assigned);
                  egap.resize( std::distance(egap.begin(), ei) );

                  /*
                  //let's resort and remove unassigned tracers
                  std::sort(egap.begin(),egap.end(),edgeg_pid_sort);

                  //remove tracers assigned to a peak
                  for(tt=0;tt<egap.size();tt++)
                    if(egap[tt].pid!=-1)
                     break;

                  //remove assigned tracers
                  egap.erase(egap.begin()+tt,egap.end());
                  */

                }//while(egap.size()>0)


                //OK, now we have to find all of the tracers
                //in the tgap array in tbuf, and update their
                //peak indices

                //first sort tbuf by index
                std::sort(tbuf.begin(),tbuf.end(),tracer_id_sort);

                //loop over tgap, and up date peak indices
                for(tt=0;tt<tgap.size();tt++)
                {
                  ia = tgap.begin() + tt;

                  //search for tgap in tbuf
                  it = std::search(tbuf.begin(),tbuf.end(),ia,ia+1,tracer_unique);

                  //once found, fix the peak_index
                  (*it).peak_index = tgap[tt].peak_index;

                }//tgap.size()

                //re-sort tbuf
                std::sort(tbuf.begin(),tbuf.end(),tracer_pid_and_density_and_id_sort);

                //we can append these shocks to 
                //the running list of blended shocks

                vector<tracer> tstmp;

                long pid = tbuf[0].peak_index;
                for(tt=0;tt<tbuf.size();tt++)
                {
                  tstmp.push_back(tbuf[tt]);


                  if(tbuf[tt].id==id_check)
                  {
                    printf("PRESENT E in tstmp pid %ld ss %ld s.id %ld tt %ld size %ld\n",tbuf[tt].id,ss,s[ss].id,tt,tbuf.size());

                    //exit(-1);
                  }

                  if(tt==tbuf.size()-1)
                  {
                    sbuf.l  = tstmp.size();
                    sbuf.o  = 0;
                    sbuf.d  = tstmp[0].d;
                    sbuf.id = tstmp[0].id;
                    set_peak_box(&sbuf,tstmp);

                    //add tracers to tmerge  //ADDTMERGE
                    smerge.push_back(sbuf);
                    for(long si=0;si<tstmp.size();si++)
                    {
                      if(tstmp[si].id==id_check)
                        printf("PRESENT TMERGE A pid %ld ss %ld s.id %ld tt %ld size %ld\n",tstmp[si].id,ss,s[ss].id,si,tstmp.size());

                      tmerge.push_back(tstmp[si]);
                    }

                    //clear tstmp
                    vector<tracer>().swap(tstmp);                    

                  }else if(tbuf[tt+1].peak_index!=pid){

                    //store this shock

                    sbuf.l  = tstmp.size();
                    sbuf.o  = 0;
                    sbuf.d  = tstmp[0].d;
                    sbuf.id = tstmp[0].id;
                    set_peak_box(&sbuf,tstmp);

                    //add tracers to tmerge //ADDTMERGE
                    smerge.push_back(sbuf);
                    for(long si=0;si<tstmp.size();si++)
                    {
                      tmerge.push_back(tstmp[si]);
                    }

                    //clear tstmp
                    vector<tracer>().swap(tstmp);

                    //update pid
                    pid = tbuf[tt+1].peak_index;
                  }

                  //make sure pid is never unassigned
                  if(pid==-1)
                  {
                    printf("ERROR pid -1\n");
                    exit(-1);
                  }
                }//loop over tbuf

                printf("************\n");
                printf("*** Blended (D) shock properties l %10ld o %10ld d %e id %10ld\n",sbuf.l,sbuf.o,sbuf.d,sbuf.id);
                printf("************\n");


                //destroy tgap and egap
                vector<tracer>().swap(tgap);
                vector<edgeg>().swap(egap);

                //free memory
                gap_tree_data.resize(extents[0][0]);
                free(gap_tree);
                flag_gap_tree = 0;

              }else{ //tgap.size()>0 -> are particles in the gap btwn peaks?

                //there are no particles in the gap

                //re-sort tbuf
                std::sort(tbuf.begin(),tbuf.end(),tracer_pid_and_density_and_id_sort);
                
                //add particles to tmerge

                vector<tracer> tstmp;

                long pid = tbuf[0].peak_index;
                for(tt=0;tt<tbuf.size();tt++)
                {
                  tstmp.push_back(tbuf[tt]);

                  if(t[s[ss].o+tt].id==id_check)
                  {
                    printf("PRESENT E in tstmp but tgap==0 ss %ld\n",ss);
                    //exit(-1);
                  }

                  if(tt==tbuf.size()-1)
                  {
                    sbuf.l  = tstmp.size();
                    sbuf.o  = 0;
                    sbuf.d  = tstmp[0].d;
                    sbuf.id = tstmp[0].id;
                    set_peak_box(&sbuf,tstmp);

                    //add tracers to tmerge  //ADDTMERGE
                    smerge.push_back(sbuf);
                    for(long si=0;si<tstmp.size();si++)
                    {
                      tmerge.push_back(tstmp[si]);
                    }

                    //clear tstmp
                    vector<tracer>().swap(tstmp);                    

                  }else if(tbuf[tt+1].peak_index!=pid){

                    //store this shock

                    sbuf.l  = tstmp.size();
                    sbuf.o  = 0;
                    sbuf.d  = tstmp[0].d;
                    sbuf.id = tstmp[0].id;
                    set_peak_box(&sbuf,tstmp);

                    //add tracers to tmerge //ADDTMERGE
                    smerge.push_back(sbuf);
                    for(long si=0;si<tstmp.size();si++)
                    {
                      tmerge.push_back(tstmp[si]);
                    }

                    //clear tstmp
                    vector<tracer>().swap(tstmp);

                    //update pid
                    pid = tbuf[tt+1].peak_index;
                  }
                }//end loop over tbuf
                printf("************\n");
                printf("*** Blended (E) shock properties l %10ld o %10ld d %e id %10ld\n",sbuf.l,sbuf.o,sbuf.d,sbuf.id);
                printf("************\n");

              }

              //Done!

//CASE 2 -- we don't find the low density peak
            }else if(it==tbuf.end()){


              //here the bounding boxes have just overlapped, and
              //the lower threshold peak does not contain the peak
              //of the higher threshold peak.
              printf("**INT #%ld: LDT DOESN'T CONTAIN HDT int peak %6ld (bs.id=%10ld) is *NOT* in ss = %10ld (s.id=%10ld (tbuf.size() %ld; bs.l %ld))\n",ti,i,(*bs)[i].id,ss,s[ss].id,tbuf.size(),(*bs)[i].l);

              //OK, what fraction of the lower density peak
              //overlaps with the higher density peak?
              long icheck = 0;

              //we can figure this out by doing a unique
              //comparison

              //add the particles in the dense peak to the 
              //a buffer to contain the union of hdt and ldt
              //tracers
              for(tt=0;tt<(*bs)[i].l;tt++)
              {
                tunion.push_back((*bt)[(*bs)[i].o+tt]);
                if((*bt)[(*bs)[i].o+tt].id==id_check)
                {
                  printf("*PRESENT BT id %ld tt %ld\n",(*bt)[(*bs)[i].o+tt].id,tt);
                  //exit(-1);
                }
              }

              for(tt=0;tt<tbuf.size();tt++)
              {
                if(tbuf[tt].id==id_check)
                {
                  printf("*PRESENT in TBUF tt %ld ss %ld s.id %ld\n",tt,ss,s[ss].id);
                  //exit(-1);
                }
              }

              //add tbuf to the end of tunion
              it = tunion.end();
              tunion.insert(it,tbuf.begin(),tbuf.end());

              //sort by id
              std::sort( tunion.begin(), tunion.end(), tracer_id_sort);

              //find the overlap from duplicated entries
              //and store them into toverlap buffer
              keep_duplicates(tunion, &toverlap);

              printf("OVERLAP size is %ld\n",toverlap.size());

             for(tt=0;tt<toverlap.size();tt++)
              {
                if(toverlap[tt].id==id_check)
                {
                  printf("*PRESENT F in OVERLAP id %ld tt %ld ss %ld s.id %ld\n",toverlap[tt].id,tt,ss,s[ss].id);
                  //exit(-1);
                }
              }
//in both union and overlap, so also in bt

              //if(toverlap.size()>0)
              if(0)
              {
                for(tt=0;tt<toverlap.size();tt++)
                  printf("overlap tt %ld id %ld d %e\n",tt,toverlap[tt].id,toverlap[tt].d);

                for(tt=0;tt<(*bs)[i].l;tt++)
                  printf("bs tt %ld id %ld d %e\n",tt,(*bt)[(*bs)[i].o+tt].id,(*bt)[(*bs)[i].o+tt].d);

              }

              //do a unique comparison on the tracer ids
              it = std::unique(tunion.begin(), tunion.end(), tracer_unique);

              //resize the union array to contain unique entries
              tunion.resize( std::distance(tunion.begin(), it) );

              printf("UNION size is %ld\n",tunion.size());


              for(tt=0;tt<tunion.size();tt++)
              {

                if(tunion[tt].id==id_check)
                {
                  printf("*PRESENT F in ORIG UNION pid %ld ss %ld s.id %ld\n",tunion[tt].id,ss,s[ss].id);
                }
              }
//still in union after restriction to unique entries

              //make a corrected list of tracers in the denser
              //peak that does not contain the overlap

              //to do that, make a list of the dense threshold peak
              //and the overlap into a buffer called tcbuf
              for(tt=0;tt<(*bs)[i].l;tt++)
                tcbuf.push_back((*bt)[(*bs)[i].o+tt]);
              
              for(tt=0;tt<toverlap.size();tt++)
                tcbuf.push_back(toverlap[tt]);

              //tcbuf contains multiple copies of
              //the duplicated particles

              printf("BEFORE SORT, tcbuf.size() %ld\n",tcbuf.size());

              //sort the ids in this corrected list
              std::sort( tcbuf.begin(), tcbuf.end(), tracer_id_sort);
              for(tt=0;tt<tcbuf.size();tt++)
              {

                if(tcbuf[tt].id==id_check)
                {
                  printf("*PRESENT F in TCBUF BEFORE TCORR pid %ld ss %ld s.id %ld\n",tcbuf[tt].id,ss,s[ss].id);
                }
              }
              //ok, we move through tcbuf and only add to tcorr those
              //elements that aren't duplicated

//duplicated as expected

              //IE we remove the duplicates from tcbuf
              keep_unique(tcbuf, &tcorr);
              printf("AFTER keep_unique, tcbuf.size() %ld tcorr.size() %ld\n",tcbuf.size(),tcorr.size());

              for(tt=0;tt<tcorr.size();tt++)
              {

                if(tcorr[tt].id==id_check)
                {
                  printf("*PRESENT F in ORIG TCORR pid %ld ss %ld s.id %ld\n",tcorr[tt].id,ss,s[ss].id);
                }
              }              

//NOT in TCORR here, since it is a duplicate

              //destroy tcbuf
              vector<tracer>().swap(tcbuf);

              //tcorr now contains the unique particles
              //in the HDT peak

              printf("*** REPEAT FOR LDT peak\n");

              //repeat for the low density threshold tracer
              for(tt=0;tt<tbuf.size();tt++)
                tcbuf.push_back(tbuf[tt]);
              
              for(tt=0;tt<toverlap.size();tt++)
                tcbuf.push_back(toverlap[tt]);

              printf("BEFORE SORT, tcbuf_ld.size() %ld\n",tcbuf.size());

              //tcbuf now contains the LDT peak, plus particles
              //duplicated with the LDT

              //sort the ids in this corrected list
              std::sort( tcbuf.begin(), tcbuf.end(), tracer_id_sort);
              for(tt=0;tt<tcbuf.size();tt++)
              {

                if(tcbuf[tt].id==id_check)
                {
                  printf("*PRESENT F in TCBUF BEFORE TCORR_LOWD pid %ld ss %ld s.id %ld\n",tcbuf[tt].id,ss,s[ss].id);
                }
              }

              //ok, we move through tcbuf and only add to tcorr_lowd those
              //elements that aren't duplicated
              keep_unique(tcbuf, &tcorr_lowd);
              printf("AFTER keep_unique, tcbuf.size() %ld tcorr_lowd.size() %ld\n",tcbuf.size(),tcorr_lowd.size());

              for(tt=0;tt<tcorr_lowd.size();tt++)
              {

                if(tcorr_lowd[tt].id==id_check)
                {
                  printf("*PRESENT F in ORIG TCORR_LOWD pid %ld ss %ld s.id %ld\n",tcorr_lowd[tt].id,ss,s[ss].id);
                }
              } 

//particle is in neither tcorr or tcorr_lowd, since it is in tbuf and in 
//bt, meaning it is in toverlap

              //build a search tree for the denser peak
              peak_tree_data.resize(extents[tcorr.size()][3]);

              //load all particle data into a tree
              for(tt=0;tt<tcorr.size();tt++)
                for(int k=0;k<3;k++)
                  peak_tree_data[tt][k] = tcorr[tt].x[k];

              //build the tree
              peak_tree = new kdtree2(peak_tree_data, true);

              //here peak_tree contains tcorr as a kdtree
              //we can search each particle in the overlap
              //to ensure it really should belong to the
              //dense peak based on the FOF search radius
              //Add to tcorr or tcorr_lowd accordingly.
              printf("BEFORE tcorr.size() %ld tcorr_lowd.size() %ld\n",tcorr.size(),tcorr_lowd.size());

              //we can only change the LDT, so
              //remove the overlap from it
              
              for(tt=0;tt<toverlap.size();tt++)
              {
                //create the query vector
                for(int k=0;k<3;k++)
                  xc[k] = toverlap[tt].x[k];

                //find any A particles within bsq
                peak_tree->r_nearest(xc,rmax,res);

                //if it belongs to the dense peak, then we'll
                //find a dense particle within the search radius
                //and if so let's add the particle from the 
                //overlap to the dense peak
                if(res.size()>0)
                {
                  if(toverlap[tt].id==id_check)
                  {
                    printf("*STORING F in ORIG TCORR pid %ld ss %ld s.id %ld\n",toverlap[tt].id,ss,s[ss].id);
                  }
                  tcorr.push_back(toverlap[tt]);
                }else{
                  //OK, this particle does not belong to the
                  //dense peak.  Add it back to the low threshold
                  //peak
                  if(toverlap[tt].id==id_check)
                  {
                    printf("*STORING F in ORIG TCORR_LOWD pid %ld ss %ld s.id %ld\n",toverlap[tt].id,ss,s[ss].id);
                  }
                  tcorr_lowd.push_back(toverlap[tt]);
                }
              }
              
              printf("AFTER  tcorr.size() %ld tcorr_lowd.size() %ld\n",tcorr.size(),tcorr_lowd.size());

              if(tcorr.size()+tcorr_lowd.size()!=tunion.size())
              //if((*bs)[i].l + tcorr_lowd.size()!=tunion.size())
              {
                //printf("ERROR LOST PARTICLE  tcorr.size() %ld tcorr_lowd.size() %ld tunion.size() %ld\n",tcorr.size(),tcorr_lowd.size(),tunion.size());
                printf("ERROR LOST PARTICLE  (*bs).l %ld tcorr_lowd.size() %ld tunion.size() %ld\n",(*bs)[i].l,tcorr_lowd.size(),tunion.size());

                fflush(stdout);
                exit(-1);
              }

              for(tt=0;tt<tcorr.size();tt++)
              {

                if(tcorr[tt].id==id_check)
                {
                  printf("*PRESENT F in ALT TCORR pid %ld peak_index %ld ss %ld s.id %ld\n",tcorr[tt].id,tcorr[tt].peak_index,ss,s[ss].id);
                }
              }   
              for(tt=0;tt<tcorr_lowd.size();tt++)
              {

                if(tcorr_lowd[tt].id==id_check)
                {
                  printf("*PRESENT F in ALT TCORR_LOWD pid %ld ss %ld s.id %ld\n",tcorr_lowd[tt].id,ss,s[ss].id);
                }
              } 

              //at this point, tcorr contains the corrected high
              //density threshold peak, adjusted to contain any
              //particles from the low density peak that may
              //have triggered the interaction

              //and tcorr_lowd contains the corrected low density
              //threshold peak, excluding any particles that
              //should properly belong to the high-density peak

              //we need to replace tbuf with tcorr_lowd

              printf("Replacing tbuf.size() %ld with tcorr_lowd.size() %ld; tcorr.size() %ld\n",tbuf.size(),tcorr_lowd.size(),tcorr.size());

              
              vector<tracer>().swap(tbuf);
              it = tbuf.end();
              tbuf.insert(it,tcorr_lowd.begin(),tcorr_lowd.end());  
              std::sort(tbuf.begin(), tbuf.end(), tracer_density_and_id_sort);

              //reset buffer information
              sbuf.l = tbuf.size();
              sbuf.o = 0;
              sbuf.d = tbuf[0].d;
              sbuf.id = tbuf[0].id;
              set_peak_box(&sbuf, tbuf);


              printf("Revised tbuf.size() %ld\n",tbuf.size());

              //destroy the peak tree
              free(peak_tree);
              peak_tree_data.resize(extents[0][0]);

              for(tt=0;tt<tunion.size();tt++)
              {

                if(tunion[tt].id==id_check)
                {
                  //printf("PRESENT F IN UNION ss %ld\n",ss);
                  printf("*PRESENT F in UNION pid %ld ss %ld s.id %ld\n",tunion[tt].id,ss,s[ss].id);

                  //exit(-1);
                }
              }

              for(tt=0;tt<tbuf.size();tt++)
              {
               if(tbuf[tt].id==id_check)
                {
                  //printf("PRESENT F IN tbuf ss %ld\n",ss);
                  printf("*PRESENT F in tbuf pid %ld ss %ld s.id %ld\n",tbuf[tt].id,ss,s[ss].id);
                }
              }

              //destroy ttest
              vector<tracer>().swap(tunion);
              vector<tracer>().swap(tcorr);
              vector<tracer>().swap(tcorr_lowd);
              vector<tracer>().swap(tcbuf);
              //vector<tracer>().swap(tout);

              //destroy toverlap, which was kept for testing
              vector<tracer>().swap(toverlap);

              //let's see if this is ever encountered.
              //exit(-1);
              flag_box_fail = 1;
              //printf("BOX FAIL %d\n",flag_box_fail);
              //exit(-1);
              printf("************\n");
              printf("*******  Blended (F) shock properties l %10ld o %10ld d %e id %10ld\n",sbuf.l,sbuf.o,sbuf.d,sbuf.id);
              printf("************\n");

              //if(toverlap.size()>0)
              //  exit(-1);



            }

            //destroy tsearch
            vector<tracer>().swap(tsearch);

          }//end loop over interactions

          //what ever is left in tbuf after
          //going through all the interactions 
          //is new, and gets added to tmerge/smerge
           
          if(tbuf.size()>0 && !flag_multi)
          {
            //sort tracers by density, then id
            std::sort(tbuf.begin(), tbuf.end(), tracer_density_and_id_sort);

            //set peak index
            for(tt=0;tt<tbuf.size();tt++)
              tbuf[tt].peak_index = tbuf[0].id;

            //reset shock properties
            sbuf.l = tbuf.size();
            sbuf.o = 0;
            sbuf.d = tbuf[0].d;
            sbuf.id = tbuf[0].id;

            //set the peak box
            set_peak_box(&sbuf, tbuf);

            //add shock and tracers to the merge lists
            smerge.push_back(sbuf);
            for(tt=0;tt<sbuf.l;tt++) //ADDTMERGE
              tmerge.push_back(tbuf[tt]);

            printf("************\n");
            printf("*** Blended (G) shock properties l %10ld o %10ld d %e id %10ld\n",sbuf.l,sbuf.o,sbuf.d,sbuf.id);
            printf("************\n");

          }
          
                    
          //destroy tbuf
          vector<tracer>().swap(tbuf);

        }//end interactions > 1


      }//end interactions>0
/*
      //check pathological case
      if(ss<s.size()-1)
      {
        if(s[ss+1].d==s[ss].d && box_collision(s[ss].min,s[ss].max,s[ss+1].min,s[ss+1].max))
        {
          printf("PATHOLOGY ss %ld id %ld d %e ss+1 %ld id %ld d %e\n",ss,s[ss].id,s[ss].d,ss+1,s[ss+1].id,s[ss+1].d);



        }

      }
*/

      //destroy the interaction list
      vector<int>().swap(interactions);

    }//end loop over low density threshold shocks (ss)
  } //end bs->size()>0

  printf("PCOUNT %ld t.size() %ld\n",pcount,t.size());


  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  //IMPORTANT BS AND BT ARE CLEARED AND REPLACED
  //////////////////////////////////////////////////


  //vacate bs and bt
  (*bs).clear();
  bt->clear();

  //bs shocks with int=1 and int>1 can
  //appear in both categories

  //tmerge + tappend + texpand can be wrong here

  //tappend is only interactions == 0
  //texpand is only interactions == 1
  //tmerge  is only interactions >  1


  printf("SIZES of tmerge %ld tappend %ld texpand %ld\n",tmerge.size(),tappend.size(),texpand.size());

  //so we need to combine into a buffer
  //and, sort by shock peak density,
  //and then store in bs and bt.
  for(ss=0;ss<smerge.size();ss++)
    sstore.push_back(smerge[ss]);
  for(tt=0;tt<tmerge.size();tt++)
  {
    if(tmerge[tt].id==id_check)
      printf("PRESENT TMERGE->TSTORE A pid %ld tt %ld size %ld\n",tmerge[tt].id,tt,tmerge.size());

    tstore.push_back(tmerge[tt]);
  }

  for(ss=0;ss<sappend.size();ss++)
    sstore.push_back(sappend[ss]);
  for(tt=0;tt<tappend.size();tt++)
    tstore.push_back(tappend[tt]);

  for(ss=0;ss<sexpand.size();ss++)
    sstore.push_back(sexpand[ss]);
  for(tt=0;tt<texpand.size();tt++)
    tstore.push_back(texpand[tt]);

  printf("BEFORE UNIQUE tstore.size() %ld t.size() %ld\n",tstore.size(),t.size());

  for(tt=0;tt<tstore.size();tt++)
  {
    if(tstore[tt].id==id_check)
      printf("PRESENT TSTORE BU A pid %ld tt %ld size %ld\n",tstore[tt].id,tt,tstore.size());
  }
  //should be able to sort by id, keep unique, then sort
  //by peak index then density then id, and reconstruct 
  //s
  //but we need to know the index of the peaks
  std::sort(tstore.begin(), tstore.end(), tracer_id_sort);

  it = std::unique(tstore.begin(), tstore.end(), tracer_unique);
  tstore.resize( std::distance(tstore.begin(), it) );

  printf("AFTER UNIQUE tstore.size() %ld t.size() %ld\n",tstore.size(),t.size());

  for(tt=0;tt<tstore.size();tt++)
  {
    if(tstore[tt].id==id_check)
      printf("PRESENT TSTORE AU A pid %ld peak_index %ld tt %ld size %ld\n",tstore[tt].id,tstore[tt].peak_index,tt,tstore.size());
  }
  //here tstore is only unique objects
  //loop, store peak_index, and then keep unique 
  //peak indices
  vector<long> pidlist;
  vector<long>::iterator pidi;
  for(tt=0;tt<tstore.size();tt++)
    pidlist.push_back(tstore[tt].peak_index);

  std::sort(pidlist.begin(),pidlist.end());
  pidi = std::unique(pidlist.begin(), pidlist.end());
  pidlist.resize( std::distance(pidlist.begin(), pidi) );

  printf("pidlist.size() %ld\n",pidlist.size());

  //pidlist contains a sorted list of unique peak indices
  //in increasing order

  //sort tracers by peak index
  std::sort(tstore.begin(), tstore.end(), tracer_pid_and_density_and_id_sort);

  vector<float> peak_densities;
  long pt = 0;
  peak_densities.push_back(tstore[0].d);
  for(tt=1;tt<tstore.size();tt++)
    if(tstore[tt].peak_index!=pidlist[pt])
    {
      peak_densities.push_back(tstore[tt].d);
      pt++;
    }

  printf("peak_densities.size() %ld\n",peak_densities.size());

  vector<tracer_key> tkey;
  tracer_key tkin;
  vector<shock>  skeep;
  vector<tracer> tkeep;

  pt = 0;
  for(tt=0;tt<tstore.size();tt++)
  {
    if(tstore[tt].peak_index!=pidlist[pt])
      pt++;
    tkin.idx = tt;
    tkin.peak_density = peak_densities[pt];
    tkin.peak_index   = tstore[tt].peak_index;
    tkin.d            = tstore[tt].d;
    tkin.id           = tstore[tt].id;

    if(tt==0)
      printf("tt %ld id %ld peak_index %ld id %ld d %e | pt %ld pid %ld pd %e\n",tt,tkin.id,tkin.peak_index,tkin.id,tkin.d,pt, pidlist[pt],peak_densities[pt]);
    tkey.push_back(tkin);
  }

  std::sort(tkey.begin(),tkey.end(), tracer_key_sort);

  //make tkeep, sorting by peak density, peak index, density, index
  for(tt=0;tt<tkey.size();tt++)
  {
    tkeep.push_back(tstore[tkey[tt].idx]);
    if(tkey[tt].id==id_check)
      printf("PRESENT TKEY A pid %ld peak_index %ld tt %ld size %ld\n",tkey[tt].id,tkey[tt].peak_index,tt,tkey.size());
  
  }

  printf("tkeep[0] id %ld peak_index %ld d %e\n",tkeep[0].id,tkeep[0].peak_index,tkeep[0].d);
  //make skeep
  long offset = 0;
  sbuf.d  = tkeep[0].d;
  sbuf.id = tkeep[0].id;
  sbuf.o  = offset;
  sbuf.l  = 1;

  for(tt=1;tt<tkeep.size();tt++)
  {
    if(tkeep[tt].peak_index!=sbuf.id)
    {
      offset += sbuf.l;
      skeep.push_back(sbuf);

      sbuf.d  = tkeep[tt].d;
      sbuf.id = tkeep[tt].id;
      sbuf.l  = 0;
      sbuf.o  = offset;
    }
    sbuf.l++;
  }
  skeep.push_back(sbuf);
  long nctot = 0;
  for(ss = 0;ss<skeep.size();ss++)
  {
    //printf("SHOCK ID LIST ss %4ld id %10ld l %6ld o %8ld d %5.4e | pid %10ld d %5.4e\n",ss,skeep[ss].id,skeep[ss].l,skeep[ss].o,skeep[ss].d,tkeep[skeep[ss].o].id,tkeep[skeep[ss].o].d);
    printf("SHOCK ID LIST ss %4ld id %8ld l %6ld o %6ld d %5.4e | pid %10ld d %5.4e\n",ss,skeep[ss].id,skeep[ss].l,skeep[ss].o,skeep[ss].d,tkeep[skeep[ss].o].id,tkeep[skeep[ss].o].d);
    nctot += skeep[ss].l;
  }
  printf("TOTAL in SHOCK LIST = %ld\n",nctot);


  printf("skeep.size() %ld tkeep.size() %ld t.size() %ld\n",skeep.size(),tkeep.size(),t.size());
  //exit(-1);

  //reset bounding boxes
  for(ss=0;ss<skeep.size();ss++)
  {
    for(int k=0;k<3;k++)
    {
      skeep[ss].min[k] =  1.0e9;
      skeep[ss].max[k] = -1.0e9;
    }

    for(tt=skeep[ss].o;tt<skeep[ss].o+skeep[ss].l;tt++)
    {
      for(int k=0;k<3;k++)
      {
        if(tkeep[tt].x[k]<skeep[ss].min[k])
          skeep[ss].min[k] = tkeep[tt].x[k];
        if(tkeep[tt].x[k]>skeep[ss].max[k])
          skeep[ss].max[k] = tkeep[tt].x[k];
      }
    }
  }

  //should be able to store and move on.

/*

            sbuf.o = 0;
            sbuf.d = tbuf[0].d;
            sbuf.id = tbuf[0].id;

  //fix offsets
  sstore[0].o = 0;
  for(ss=1;ss<sstore.size();ss++)
    sstore[ss].o = sstore[ss-1].o + sstore[ss-1].l;

  //there might be duplicates, because of marginal 
  //overlapping shocks

  //sort shocks by peak id
  std::sort(sstore.begin(), sstore.end(), shock_id_sort);

*/



  //check for duplicates

//  *****************************

  //create tracer keys for sorting
/*
  vector<long> pidlist;
  vector<long>::iterator pidi;
  for(tt=0;tt<tstore.size();tt++)
    pidlist.push_back(tstore[tt].peak_index);


  vector<tracer_key> tkey;
  vector<float>  skey;

  shock_key skin;
  for(ss=0;ss<sstore.size())

  tracer_key tkin;
  tt = 0;
  for(ss=0;ss<sstore.size();ss++)
  {
    for(tt=sstore[ss].o;tt<sstore[ss].osstore[ss].l;tidx++)
    tkin.idx = tt;
    tkin.peak_density = tstore[tstore[tt].peak_idx].

  }

  //let's just make sure tstore reduces to t
  //keep unique entries
  std::sort(tstore.begin(), tstore.end(), tracer_id_sort);
  it = std::unique(tstore.begin(), tstore.end(), tracer_unique);
  tstore.resize( std::distance(tstore.begin(), it) );


  //sort by peak, then density, then id
  std::sort(tstore.begin(), tstore.end(), tracer_pid_and_density_and_id_sort);
*/
  //for(ss=0;ss<store.size();ss++)
    //printf("SORTING SCHEME SS ss %ld id %ld d %e\n",ss,store[ss].id,store[ss].d);
  //HEREHEREHERE

/*
  tt=0;
  ss=0;
  long llcheck=1;
  printf("SORTING SCHEME TT tt %ld id %ld pid %ld d %e ss %ld l %ld o %ld llcheck %ld id %ld d %e\n",tt,tstore[tt].id,tstore[tt].peak_index,tstore[tt].d,ss,sstore[ss].l,sstore[ss].o,llcheck,sstore[ss].id,sstore[ss].d);
  for(tt=1;tt<tstore.size();tt++)
  {
    if(tstore[tt].peak_index!=tstore[tt-1].peak_index)
    {
      ss++;
      printf("SORTING SCHEME TT tt %ld id %ld pid %ld d %e ss %ld l %ld o %ld lcheck %ld id %ld d %e\n",tt,tstore[tt].id,tstore[tt].peak_index,tstore[tt].d,ss,sstore[ss].l,sstore[ss].o,llcheck,sstore[ss].id,sstore[ss].d);
      llcheck = 0;
    }else{
      llcheck++;
    }

  }

//exit(-1);
//  *****************************

  //first, keep the first peak
  skeep.push_back(sstore[0]);
  for(tt=sstore[0].o;tt<sstore[0].o+sstore[0].l;tt++)
    tkeep.push_back(tstore[tt]);

  //now loop over peaks
  for(ss=1;ss<sstore.size();ss++)
  {
    //if the peak is not a duplicate,
    //keep it
    if(sstore[ss].id!=sstore[ss-1].id)
    {
      skeep.push_back(sstore[ss]);
    }else{
      skeep[skeep.size()-1].l += sstore[ss].l;
    }

    //but be sure to store it's particles
    //in any case
    for(tt=sstore[ss].o;tt<sstore[ss].o+sstore[ss].l;tt++)
      tkeep.push_back(tstore[tt]);
  }

  //tkeep / tstore can be larger than t here
  //but tkeep recovers all particles in tstore
  printf("SIZES t %ld tkeep %ld tstore %ld\n",t.size(),tkeep.size(),tstore.size());
  //cin.get();


  printf("SECOND SIZES t %ld tkeep %ld tstore %ld\n",t.size(),tkeep.size(),tstore.size());
  //cin.get();

  //indeed tstore usually reduces to t here under unique retainment
  //but when there is a problem it doesn't.
  //PCOUNT 722488 t.size() 722488
  //SIZES t 722488 tkeep 746573 tstore 746573
  //SECOND SIZES t 722488 tkeep 746573 tstore 722452
  //****** 13 TBUF before unique 622 after 380 Delta 242 skeep.size() 622
  //****** 64 TBUF before unique 1244 after 735 Delta 509 skeep.size() 1244
  //****** 116 TBUF before unique 6413 after 3582 Delta 2831 skeep.size() 6413 
  //****** 143 TBUF before unique 18121 after 10579 Delta 7542 skeep.size() 18121
  //****** 208 TBUF before unique 15691 after 8420 Delta 7271 skeep.size() 15691
  //****** 264 TBUF before unique 11001 after 6232 Delta 4769 skeep.size() 11001
  //****** 590 TBUF before unique 890 after 501 Delta 389 skeep.size() 890
  //****** 593 TBUF before unique 841 after 397 Delta 444 skeep.size() 841
  //****** 808 TBUF before unique 329 after 205 Delta 124 skeep.size() 329
  //LENGTH ERROR LOST TRACERS t.size 722488 tfinal.size() 722452


  //reset the offsets
  skeep[0].o = 0;
  for(ss=1;ss<skeep.size();ss++)
  {
    skeep[ss].o = skeep[ss-1].o + skeep[ss-1].l;
  }

  //need to remove duplicate tracers
  vector<tracer> tfinal, tfbuf;

  //loop over the shocks that we're keeping
  for(ss=0;ss<skeep.size();ss++)
  {
    //correct skeep
    printf("SKEEP ss %ld l %6ld o %6ld id %10ld d %e\n",ss,skeep[ss].l,skeep[ss].o,skeep[ss].id,skeep[ss].d);

    //first, put the tracers from
    //this shock into a buffer
    for(tt=skeep[ss].o;tt<skeep[ss].o+skeep[ss].l;tt++)
      tfbuf.push_back(tkeep[tt]);

    //printf("BEFORE T SORT BY ID tbuf[0].id %ld\n",tfbuf[0].id);

    //sort by id
    std::sort(tfbuf.begin(),tfbuf.end(),tracer_id_sort);

    //keep unique entries
    tbb = tfbuf.size();
    it = std::unique(tfbuf.begin(), tfbuf.end(), tracer_unique);
    tfbuf.resize( std::distance(tfbuf.begin(), it) );

    //resort tfbuf by density
    std::sort(tfbuf.begin(),tfbuf.end(),tracer_density_and_id_sort);
    tba = tfbuf.size();

    if(tba!=tbb)
      printf("****** %ld TBUF before unique %ld after %ld Delta %ld skeep.size() %ld\n",ss,tbb,tba,tbb-tba,skeep[ss].l);

    //put tbuf in tfinal
    for(tt=0;tt<tfbuf.size();tt++)
    { 
      if(tfbuf[tt].id==1626042 || tfbuf[tt].peak_index==1467809)
        printf("TBUF INTO FINAL tt %ld id %ld peak_index %ld d %e\n",tt,tfbuf[tt].id,tfbuf[tt].peak_index,tfbuf[tt].d);
      tfinal.push_back(tfbuf[tt]);
    }

    //correct skeep length
    skeep[ss].l = tfbuf.size();

    //check for an error
    if(skeep[ss].id!=tfbuf[0].id)
    {
      printf("ID ERROR ss %ld skeep.id %ld tfbuf.id %ld skeep.l %ld\n",ss,skeep[ss].id,tfbuf[0].id,skeep[ss].l);
      exit(-1);
    }

    //destroy tfbuf
    vector<tracer>().swap(tfbuf);
  }

  //correct offsets
  skeep[0].o = 0;
  for(ss=1;ss<skeep.size();ss++)
  {
    skeep[ss].o = skeep[ss-1].o + skeep[ss-1].l;
  }

  //OK, now sort skeep based on density
  std::sort(skeep.begin(), skeep.end(), shock_density_sort);

  //there is a pathological case where
//#error deal with case where a proximate equal density peak appears in int=1, 
//#error but you've already enshrined this peak and so it's duplicated in 
//#error the interactions==1 case previously.

  //what does skeep look like now?
  for(ss=0;ss<skeep.size();ss++)
    printf("MAKING TFINAL ss %ld skeep.id %ld skeep.o %ld skeep.e %e\n",ss,skeep[ss].id,skeep[ss].o,skeep[ss].d);
*/

  //Store the tracers in blended shocks
  /*
  tt=0;
  for(ss=0;ss<skeep.size();ss++)
    for(tt=skeep[ss].o;tt<skeep[ss].o+skeep[ss].l;tt++)
      bt->push_back(tkeep[tt]);
      //bt->push_back(tfinal[tt]);
*/


  //Store the tracers in blended shocks
  for(tt=0;tt<tkeep.size();tt++)
    bt->push_back(tkeep[tt]);
  

/*
  //fix offsets for the last time, after sorting tracers
  skeep[0].o = 0;
  for(ss=1;ss<skeep.size();ss++)
  {
    skeep[ss].o = skeep[ss-1].o + skeep[ss-1].l;
  }
*/
  
  //put in shocks
  for(ss=0;ss<skeep.size();ss++)
    bs->push_back(skeep[ss]);

  //printf("CHECKING KEEP *********\n");
  //for(ss=0;ss<skeep.size();ss++)
    //printf("ss %6ld KEEP l %6ld o %6ld id %10ld d %e\n",ss,skeep[ss].l,skeep[ss].o,skeep[ss].id,skeep[ss].d);

  //OK, now skeep and tkeep are our 
  //corrected shocks. But tkeep probably
  //contains duplicates 

  //sort sstore based on density
  //std::sort(sstore.begin(), sstore.end(), shock_density_sort);
/*

  vector<shock>  skeep;
  vector<tracer> tkeep;
  //check for duplicates and remove them
  for(ss=0;ss<sstore.size();ss++)
  {
    if(ss==sstore.size()-1)
    {
      skeep.push_back(sstore[ss]);
      for(tt=sstore[ss].o;tt<sstore[ss].l;tt++)
        tkeep.push_back(tstore[tt]);
    }else if(sstore[ss].id==sstore[ss+1].id)
    {
      //we have a duplicate
      skeep.push_back(sstore[ss]);
      for(tt=sstore[ss].o;tt<sstore[ss].l;tt++)
        tkeep.push_back(tstore[tt]);

      ss++; //skip duplicate
    }else{
      skeep.push_back(sstore[ss]);
      for(tt=sstore[ss].o;tt<sstore[ss].l;tt++)
        tkeep.push_back(tstore[tt]);
    }
  }

  //
  //fix offsets
  skeep[0].o = 0;
  for(ss=1;ss<skeep.size();ss++)
    skeep[ss].o = skeep[ss-1].o + skeep[ss-1].l;

  //put tracers in bt according to new order
  tt=0;
  for(ss=0;ss<skeep.size();ss++)
    for(tt=skeep[ss].o;tt<skeep[ss].o+skeep[ss].l;tt++)
      bt->push_back(tkeep[tt]);
    */
  
  /*
  tt=0;
  for(ss=0;ss<sstore.size();ss++)
    for(tt=sstore[ss].o;tt<sstore[ss].o+sstore[ss].l;tt++)
      bt->push_back(tstore[tt]);
  


  //fix offsets for the last time, after sorting tracers
  sstore[0].o = 0;
  for(ss=1;ss<sstore.size();ss++)
  {
    sstore[ss].o = sstore[ss-1].o + sstore[ss-1].l;
  }
  

  //put in shocks
  for(ss=0;ss<sstore.size();ss++)
    bs->push_back(sstore[ss]);
  */

//  for(ss=0;ss<skeep.size();ss++)
  //  bs->push_back(skeep[ss]);



  //print information about the shocks.

  shock sbufB;
  printf("********** Blended Shock List\n");
  for(ss=0;ss<(*bs).size();ss++)
  {
    sbuf = (*bs)[ss];


      
    if(ss<scomp.size())
    {
      sbufB = scomp[ss];
      printf("ss %6ld l %6ld o %6ld d %e id %10ld\t\tPREV ss %6ld l %6ld o %6ld d %e id %ld\n",ss,sbuf.l,sbuf.o,sbuf.d,sbuf.id,ss,sbufB.l,sbufB.o,sbufB.d,sbufB.id);
    }else{
      printf("ss %6ld l %6ld o %6ld d %e id %10ld\n",ss,sbuf.l,sbuf.o,sbuf.d,sbuf.id);
    }

    if(ss>0)
      if(sbuf.id==(*bs)[ss-1].id)
      {
        printf("DUPLICATE SHOCKS ERROR\n");
        exit(-1);
      }
  }
  printf("********** END Blended Shock List\n");



  //if necessary, free memory
  //associated with tree
  if(flag_tree_build)
  {
    bs_data.resize(extents[0][0]);
    free(bs_tree);
  }

/*
  //for debugging
  //seems ok
  //if(flag_box_fail)
    //exit(-1);
  for(tt=0;tt<tfinal.size();tt++)
    if((tfinal[tt].id==294744)||(tfinal[tt].id==16470779)||(tfinal[tt].id==43349505)||(tfinal[tt].id==22663579))
     printf("tt   %ld t id %ld peak_index %ld d %e\n",tt,tfinal[tt].id,tfinal[tt].peak_index,tfinal[tt].d);

*/


  //LENGTH CHECK
  //all the particles in t should be in tfinal
  //if(t.size()!=tfinal.size())
  if(t.size()!=tkeep.size())
  {
    //printf("LENGTH ERROR LOST TRACERS t.size %ld tfinal.size() %ld\n",t.size(),tfinal.size());
    printf("LENGTH ERROR LOST TRACERS t.size %ld tfinal.size() %ld\n",t.size(),tkeep.size());

/*
    for(tt=0;tt<tfinal.size();tt++)
    {
      if(tt<tfinal.size()-1)
        if(tfinal[tt].peak_index!=1467809)
          if(tfinal[tt+1].peak_index == 1467809)
            printf("BEFORE LIST tt   %ld t id %ld peak_index %ld d %e\n",tt,tfinal[tt].id,tfinal[tt].peak_index,tfinal[tt].d);

      if(tfinal[tt].peak_index == 1467809)
        printf("LIST tt   %ld t id %ld peak_index %ld d %e\n",tt,tfinal[tt].id,tfinal[tt].peak_index,tfinal[tt].d);
    }
    */


    //OK, we have a problem.  Let's sort both t and tfinal by tracer ID
    //then we can search and identify what is missing.
    std::sort(t.begin(), t.end(), tracer_id_sort);
    //std::sort(tfinal.begin(),tfinal.end(),tracer_id_sort);

    std::sort(tkeep.begin(),tkeep.end(),tracer_id_sort);

    for(ss=0;ss<s.size();ss++)
    {
      if((s[ss].id==16470779)||(s[ss].id==22663579)||(s[ss].id==65341979))
      {
        printf("ss %ld s.id %ld s.l %ld s.d %e\n",ss,s[ss].id,s[ss].l,s[ss].d);
      }
    }

    for(tt=0;tt<t.size()-1;tt++)
    {
      //if(t[tt].id!=tfinal[tt].id)
      if(t[tt].id!=tkeep[tt].id)
      {
        /*
        printf("tt-2 %ld t id %ld peak_index %ld tfinal %ld peak_index %ld\n",tt-2,t[tt-2].id,t[tt-2].peak_index,tfinal[tt-2].id,tfinal[tt-2].peak_index);
        printf("tt-1 %ld t id %ld peak_index %ld tfinal %ld peak_index %ld\n",tt-1,t[tt-1].id,t[tt-1].peak_index,tfinal[tt-1].id,tfinal[tt-1].peak_index);
        printf("tt   %ld t id %ld peak_index %ld tfinal %ld peak_index %ld\n",tt,t[tt].id,t[tt].peak_index,tfinal[tt].id,tfinal[tt].peak_index);
        printf("tt+1 %ld t id %ld peak_index %ld tfinal %ld peak_index %ld\n",tt+1,t[tt+1].id,t[tt+1].peak_index,tfinal[tt+1].id,tfinal[tt+1].peak_index);         
        printf("tt+2 %ld t id %ld peak_index %ld tfinal %ld peak_index %ld\n",tt+2,t[tt+2].id,t[tt+2].peak_index,tfinal[tt+2].id,tfinal[tt+2].peak_index);         
        printf("tt+3 %ld t id %ld peak_index %ld tfinal %ld peak_index %ld\n",tt+3,t[tt+3].id,t[tt+3].peak_index,tfinal[tt+3].id,tfinal[tt+3].peak_index);
        printf("tt+4 %ld t id %ld peak_index %ld tfinal %ld peak_index %ld\n",tt+4,t[tt+4].id,t[tt+4].peak_index,tfinal[tt+4].id,tfinal[tt+4].peak_index);
        */
        printf("tt-2 %ld t id %ld peak_index %ld tfinal %ld peak_index %ld\n",tt-2,t[tt-2].id,t[tt-2].peak_index,tkeep[tt-2].id,tkeep[tt-2].peak_index);
        printf("tt-1 %ld t id %ld peak_index %ld tfinal %ld peak_index %ld\n",tt-1,t[tt-1].id,t[tt-1].peak_index,tkeep[tt-1].id,tkeep[tt-1].peak_index);
        printf("tt   %ld t id %ld peak_index %ld tfinal %ld peak_index %ld\n",tt,t[tt].id,t[tt].peak_index,tkeep[tt].id,tkeep[tt].peak_index);
        printf("tt+1 %ld t id %ld peak_index %ld tfinal %ld peak_index %ld\n",tt+1,t[tt+1].id,t[tt+1].peak_index,tkeep[tt+1].id,tkeep[tt+1].peak_index);         
        printf("tt+2 %ld t id %ld peak_index %ld tfinal %ld peak_index %ld\n",tt+2,t[tt+2].id,t[tt+2].peak_index,tkeep[tt+2].id,tkeep[tt+2].peak_index);         
        printf("tt+3 %ld t id %ld peak_index %ld tfinal %ld peak_index %ld\n",tt+3,t[tt+3].id,t[tt+3].peak_index,tkeep[tt+3].id,tkeep[tt+3].peak_index);
        printf("tt+4 %ld t id %ld peak_index %ld tfinal %ld peak_index %ld\n",tt+4,t[tt+4].id,t[tt+4].peak_index,tkeep[tt+4].id,tkeep[tt+4].peak_index);        
        break;
      }

    }
    //5708317 is missing

    exit(-1);
  }else{
    //printf("PASS ALL TRACERS FOUND t.size %ld tfinal.size() %ld\n",t.size(),tfinal.size());
    printf("PASS ALL TRACERS FOUND t.size %ld tfinal.size() %ld\n",t.size(),tkeep.size());

  }

  //destroy buffer memory
  vector<tracer>().swap(tmerge);
  vector<tracer>().swap(tappend);
  vector<tracer>().swap(texpand);

  vector<tracer>().swap(tstore);
  vector<tracer>().swap(tkeep);
  //vector<tracer>().swap(tfinal);

  vector<shock>().swap(smerge);
  vector<shock>().swap(sappend);
  vector<shock>().swap(sexpand);
  vector<shock>().swap(sstore);
  vector<shock>().swap(skeep);

  //destroy shock comparison list
  vector<shock>().swap(scomp);
}
