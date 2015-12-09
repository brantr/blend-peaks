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

  vector<shock>  sappend; //new shocks to add to the list
  vector<tracer> tappend;

  vector<shock>  sexpand; //shocks that need to be expanded
  vector<tracer> texpand;

  vector<shock>  smerge; //shocks that need to be merged
  vector<tracer> tmerge;

  shock  sbuf;            //buffer shocks
  vector<tracer> tbuf;    //buffer tracers
  vector<tracer> tsearch; //buffer tracers


  vector<tracer> tunion;    //union of high and low rho threshold peaks
  vector<tracer> toverlap;  //overlap between high and low rho threshold peaks
  vector<tracer> tcorr;       //corrected high-density threshold peak
  vector<tracer> tcorr_lowd;  //corrected low density threshold peak
  vector<tracer> tcbuf;


  vector<int> interactions;
  vector<int> n_append;
  vector<int> n_merge;
  vector<int> n_expand;

  vector<tracer>::iterator it;
  vector<tracer>::iterator ia;

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

    //loop over shocks 
    for(ss=0;ss<s.size();ss++)
    {
      //add all the tracers from
      //this shock to a buffer
      for(tt=0;tt<s[ss].l;tt++)
        tbuf.push_back(t[s[ss].o+tt]);
      
      //sort all the tracers by density, and then by id
      std::sort(tbuf.begin(), tbuf.end(), tracer_density_and_id_sort);

      //add to blended tracers
      for(tt=0;tt<tbuf.size();tt++)
      {
        //set peak index
        tbuf[tt].peak_index = tbuf[0].id;

        //add to blended tracers
        bt->push_back(tbuf[tt]);
      }

      //reset peak index
      s[ss].id = tbuf[0].id;

      //set the peak box
      set_peak_box(&s[ss], tbuf);

      //add peak to peak list
      bs->push_back(s[ss]);

      //destroy tbuf
      vector<tracer>().swap(tbuf);

    }//end loop over shocks

  //end bs->size()==0
  }else{

    //this is our main work loop
    //and bs->size()>0

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

        //Put tracers into a buffer
        for(tt=0;tt<s[ss].l;tt++)
          tbuf.push_back(t[s[ss].o+tt]);
      
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

        printf("Added ss %ld (id = %ld) to the append list.\n",ss,s[ss].id);

        //destroy tracer buffer
        vector<tracer>().swap(tbuf);

      }//end interactions==0

      //ok, there is at least one interaction
      if(interactions.size()>0)
      {

        //if there is only one interaction, the shock
        //just needs to be expanded.  Let's get some info

        //note for *one* interaction, the peak with the lower
        //density threshold has to be the same as the previously
        //identified peak.  We don't need to be any more careful.
        //note that tmerge and smerge contain the lists of the
        //simply grown shocks
        if(interactions.size()==1)
        {
          i = interactions[0];
          printf("Blended shock box %e %e %e %e %e %e\n",(*bs)[i].min[0],(*bs)[i].min[1],(*bs)[i].min[2],(*bs)[i].max[0],(*bs)[i].max[1],(*bs)[i].max[2]);
          printf("New shock box     %e %e %e %e %e %e\n",s[ss].min[0],s[ss].min[1],s[ss].min[2],s[ss].max[0],s[ss].max[1],s[ss].max[2]);
          printf("Blended shock properties l %10ld o %10ld d %e id %10ld\n",(*bs)[i].l,(*bs)[i].o,(*bs)[i].d,(*bs)[i].id);
          printf("New shock properties     l %10ld o %10ld d %e id %10ld\n",s[ss].l,s[ss].o,s[ss].d,s[ss].id);

          //add the high and low density threshold
          //shocks' tracers to a buffer
          for(tt=0;tt<s[ss].l;tt++)
            tbuf.push_back(t[s[ss].o+tt]);
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
          sbuf.o  = tmerge.size();
          sbuf.d  = tbuf[0].d;
          sbuf.id = tbuf[0].id;


          //set the peak box
          set_peak_box(&sbuf, tbuf);

          //adjust the tracers' peak indices
          for(tt=0;tt<tbuf.size();tt++)
            tbuf[tt].peak_index = sbuf.id;

          //append tbuf to tmerge
          it = tmerge.end();
          tmerge.insert(it,tbuf.begin(),tbuf.end());

          //append sbuf to smerge
          smerge.push_back(sbuf);

          //destroy tbuf
          vector<tracer>().swap(tbuf);

        }//end interactions==1


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
          if(!flag_tree_build) //build bs_tree
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
          }

          //sort by density, then id
          std::sort(tbuf.begin(), tbuf.end(), tracer_density_and_id_sort);

          //reset group id
          s[ss].id = tbuf[0].id;

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


            //Here, we need to consider separately
            //whether the interactions are merges,
            //or if the bounding boxes have just overlapped
            //needs to be fixed!!!!!!
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
                //error for same, check id; 
                //otherwise build second tree and then adjust peak_index to assign

                /*
                for(long k=0;k<res.size();k++)
                  printf("ss %ld tt %ld k %ld res.dis() %e\n",ss,tt,k,res[k].dis);
                */
              }
            }else{
              //here the bounding boxes have just overlapped, and
              //the lower threshold peak does not contain the peak
              //of the higher threshold peak.
              printf("interaction peak %6ld (id=%10ld) is *NOT* in ss = %10ld (id=%10ld)\n",i,(*bs)[i].id,ss,s[ss].id);
              //exit(-1);

              //OK, what fraction of the lower density peak
              //overlaps with the higher density peak?
              long icheck = 0;

              //we can figure this out by doing a unique
              //comparison

              //add the particles in the dense peak to the 
              //a buffer to contain the union of hdt and ldt
              //tracers
              for(tt=0;tt<(*bs)[i].l;tt++)
                tunion.push_back((*bt)[(*bs)[i].o+tt]);

              //append tbuf to the end of tunion
              it = tunion.end();
              tunion.insert(it,tbuf.begin(),tbuf.end());

              //sort by id
              std::sort( tunion.begin(), tunion.end(), tracer_id_sort);

              //find the overlap from duplicated entries
              //and store them into toverlap buffer
              keep_duplicates(tunion, &toverlap);

              printf("OVERLAP size is %ld\n",toverlap.size());

              //do a unique comparison on the tracer ids
              it = std::unique(tunion.begin(), tunion.end(), tracer_unique);

              //resize the union array to contain unique entries
              tunion.resize( std::distance(tunion.begin(), it) );

              printf("UNION size is %ld\n",tunion.size());


/*
              //////////////////////////////
              //output for checking

              //save the original denser peak
              vector<tracer> tout;
              for(tt=(*bs)[i].o;tt<(*bs)[i].o+(*bs)[i].l;tt++)
                tout.push_back((*bt)[tt]);
              std::sort( tout.begin(), tout.end(), tracer_position);

              FILE *fpcheck;
              fpcheck = fopen("test_bsl.txt","w");
              for(tt=0;tt<tout.size();tt++)
                fprintf(fpcheck,"%e\t%e\t%e\t%ld\n",tout[tt].x[0],tout[tt].x[1],tout[tt].x[2],tout[tt].id);
              fclose(fpcheck);
              

              // allparticles in low density peak
              std::sort( tbuf.begin(), tbuf.end(), tracer_position);

              fpcheck = fopen("test_tbuf.txt","w");
              //fprintf(fpcheck,"%ld\n",tbuf.size());
              for(tt=0;tt<tbuf.size();tt++)
                fprintf(fpcheck,"%e\t%e\t%e\t%ld\n",tbuf[tt].x[0],tbuf[tt].x[1],tbuf[tt].x[2],tbuf[tt].id);;
              fclose(fpcheck);

              // all particles across both peaks
              std::sort( tunion.begin(), tunion.end(), tracer_position);

              fpcheck = fopen("test_tunion.txt","w");
              for(tt=0;tt<tunion.size();tt++)
                fprintf(fpcheck,"%e\t%e\t%e\t%ld\n",tunion[tt].x[0],tunion[tt].x[1],tunion[tt].x[2],tunion[tt].id);
              fclose(fpcheck);


              // all particles in in overlap

              std::sort( toverlap.begin(), toverlap.end(), tracer_position);
              fpcheck = fopen("test_toverlap.txt","w");
              for(tt=0;tt<toverlap.size();tt++)
                fprintf(fpcheck,"%e\t%e\t%e\t%ld\n",toverlap[tt].x[0],toverlap[tt].x[1],toverlap[tt].x[2],toverlap[tt].id);
              fclose(fpcheck);

              //END output for checking
              //////////////////////////////
*/

              //make a corrected list of tracers in the denser
              //peak that does not contain the overlap

              //to do that, make a list of the dense threshold peak
              //and the overlap into a buffer called tcbuf
              for(tt=0;tt<(*bs)[i].l;tt++)
                tcbuf.push_back((*bt)[(*bs)[i].o+tt]);
              
              for(tt=0;tt<toverlap.size();tt++)
                tcbuf.push_back(toverlap[tt]);

              //sort the ids in this corrected list
              std::sort( tcbuf.begin(), tcbuf.end(), tracer_id_sort);

              //ok, we move through tcbuf and only add to tcorr those
              //elements that aren't duplicated
              keep_unique(tcbuf, &tcorr);
/*
              tt=0;
              while(tt<tcbuf.size()-1)
              {
                if(tcbuf[tt].id!=tcbuf[tt+1].id)
                {
                  tcorr.push_back(tcbuf[tt]);
                }else{
                  tt++;
                }
                tt++;
              }
              if(tcbuf[tcbuf.size()-2].id!=tcbuf[tcbuf.size()-1].id)
                tcorr.push_back(tcbuf[tcbuf.size()-1]);
*/
              //destroy tcbuf
              vector<tracer>().swap(tcbuf);

              //repeat for the low density threshold tracer
              for(tt=0;tt<tbuf.size();tt++)
                tcbuf.push_back(tbuf[tt]);
              
              for(tt=0;tt<toverlap.size();tt++)
                tcbuf.push_back(toverlap[tt]);

              //sort the ids in this corrected list
              std::sort( tcbuf.begin(), tcbuf.end(), tracer_id_sort);

              //ok, we move through tcbuf and only add to tcorr those
              //elements that aren't duplicated
              tt=0;
              while(tt<tcbuf.size()-1)
              {
                if(tcbuf[tt].id!=tcbuf[tt+1].id)
                {
                  tcorr_lowd.push_back(tcbuf[tt]);
                }else{
                  tt++;
                }
                tt++;
              }
              if(tcbuf[tcbuf.size()-2].id!=tcbuf[tcbuf.size()-1].id)
                tcorr_lowd.push_back(tcbuf[tcbuf.size()-1]);




              /*
              for(tt=0;tt<tcorr.size();tt++)
              {
                printf("TCORR AFTER SORT tt %ld id %ld\n",tt,tcorr[tt].id);
              }

              fpcheck = fopen("test_tcorr.after_sort.txt","w");
              //fprintf(fpcheck,"%ld\n",ttest.size());
              for(tt=0;tt<tcorr.size();tt++)
                fprintf(fpcheck,"%e\t%e\t%e\t%ld\n",tcorr[tt].x[0],tcorr[tt].x[1],tcorr[tt].x[2],tcorr[tt].id);
              fclose(fpcheck);

              //do a unique comparison on the tracer ids
              //it = std::unique(tcorr.begin(), tcorr.end(), tracer_unique);

              //resize the search array
              //tcorr.resize( std::distance(tcorr.begin(), it) ); 
              /*
              ia = std::adjacent_find(tsearch.begin(), tsearch.end(), tracer_unique);
              if(ia!=tsearch.end())
              {
                ttest.push_back(*ia);
              
                while(ia!=tsearch.end())
                {
                  ia = std::adjacent_find(++ia, tsearch.end(), tracer_unique);
                  if(ia!=tsearch.end())
                    ttest.push_back(*ia);
                }
              }
              */     

/*
              //////////////////////////////////
              //BEGIN output for checking

              //trick sorting
              std::sort( tcorr.begin(), tcorr.end(), tracer_position);

              fpcheck = fopen("test_tcorr.txt","w");
              //fprintf(fpcheck,"%ld\n",ttest.size());
              for(tt=0;tt<tcorr.size();tt++)
                fprintf(fpcheck,"%e\t%e\t%e\t%ld\n",tcorr[tt].x[0],tcorr[tt].x[1],tcorr[tt].x[2],tcorr[tt].id);
              fclose(fpcheck);

              //END output for checking
              //////////////////////////////////
*/

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
              printf("BEFORE tcorr.size() %ld tcorr_lowd.size() %ld\n",tcorr.size(),tcorr_lowd.size());
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
                  tcorr.push_back(toverlap[tt]);
                }else{
                  //OK, this particle does not belong to the
                  //dense peak.  Add it back to the low threshold
                  //peak
                  tcorr_lowd.push_back(toverlap[tt]);
                }
              }
              printf("AFTER  tcorr.size() %ld tcorr_lowd.size() %ld\n",tcorr.size(),tcorr_lowd.size());

              if(tcorr.size()+tcorr_lowd.size()!=tunion.size())
              {
                printf("ERROR LOST PARTICLE  tcorr.size() %ld tcorr_lowd.size() %ld tunion.size() %ld\n",tcorr.size(),tcorr_lowd.size(),tunion.size());
                fflush(stdout);
                exit(-1);
              }

              exit(-1);

              //at this point, tcorr contains the corrected high
              //density threshold peak, adjusted to contain any
              //particles from the low density peak that may
              //have triggered the interaction

              //and tcorr_lowd contains the corrected low density
              //threshold peak, excluding any particles that
              //should properly belong to the high-density peak

//#error   DECIDE HOW TO INCORPORATE THESE INTO CATALOGUES


              //destroy the peak tree
              free(peak_tree);
              peak_tree_data.resize(extents[0][0]);

              //destroy ttest
              vector<tracer>().swap(tunion);
              vector<tracer>().swap(toverlap);
              vector<tracer>().swap(tunion);
              vector<tracer>().swap(tcorr);
              vector<tracer>().swap(tcorr_lowd);
              vector<tracer>().swap(tcbuf);
              //vector<tracer>().swap(tout);

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

        }//end interactions > 1

      }//end interactions>0

      //destroy the interaction list
      vector<int>().swap(interactions);

    }//end loop over low density threshold shocks (ss)

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

  } //end bs->size()>0

  //if necessary, free memory
  //associated with tree
  if(flag_tree_build)
  {
    bs_data.resize(extents[0][0]);
    free(bs_tree);
  }
}