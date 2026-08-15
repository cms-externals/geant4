
#include "../test48/PlotMuMinusStopping.C"

// void MuMinusRegre( char target[3]="Al" )
void MuMinusRegre( std::string tg="Al" )
{

   // std::string tg(target);
   
   for ( int m=0; m<NModelsMuons; ++m )
   {
      PlotMuMinusRegre( tg, ModelsMuons[m] );
   }
      
   return;

}
