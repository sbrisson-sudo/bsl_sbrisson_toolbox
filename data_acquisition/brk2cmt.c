/*
 Program that translates the unique 8 character alphanumeric 
 identifier to the original 14 character CMT name.
 V.Lekic June 1, 2006
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#define PIE acos(-1.)


// K.Yuan Sept. 12, 2014, main -> int main
int main(int argc,char *argv[]) {

 char cmtname[15],brkname[9];
 int century, year1, year2,  month, day, hour, minute, helper, faktor, remaindr;
 char datatype, lastletter, helper2[2], yrstrng[2],monthnm[2],minutnm[2];

 if(argc!=2)
// K.Yuan Sept. 12, 2014, usage() -> print ...
 {
   printf("usage: brk2cmt brkcode(a8)\n");
   return 0;
 }
// usage("brk2cmt brkcode(a8)");    
 strncpy(brkname,argv[1],8);
 brkname[8]='\0';

 if(brkname[0]=='A') {cmtname[0]='B'; cmtname[1]='1'; cmtname[2]='9';}
 if(brkname[0]=='B') {cmtname[0]='B'; cmtname[1]='2'; cmtname[2]='0';}
 if(brkname[0]=='C') {cmtname[0]='C'; cmtname[1]='1'; cmtname[2]='9';}
 if(brkname[0]=='D') {cmtname[0]='C'; cmtname[1]='2'; cmtname[2]='0';}
 if(brkname[0]=='E') {cmtname[0]='M'; cmtname[1]='1'; cmtname[2]='9';}
 if(brkname[0]=='F') {cmtname[0]='M'; cmtname[1]='2'; cmtname[2]='0';}
 if(brkname[0]=='G') {cmtname[0]='S'; cmtname[1]='1'; cmtname[2]='9';}
 if(brkname[0]=='H') {cmtname[0]='S'; cmtname[1]='2'; cmtname[2]='0';}
 cmtname[3]=brkname[1];
 cmtname[4]=brkname[2];

 if(brkname[5]=='0') faktor=0; 
 if(brkname[5]=='1') faktor=1; 
 if(brkname[5]=='2') faktor=2; 
 if(brkname[5]=='3') faktor=3; 
 if(brkname[5]=='4') faktor=4; 
 if(brkname[5]=='5') faktor=5; 
 if(brkname[5]=='6') faktor=6; 
 if(brkname[5]=='7') faktor=7; 
 if(brkname[5]=='8') faktor=8; 
 if(brkname[5]=='9') faktor=9; 
 if(brkname[5]=='A') faktor=10; 
 if(brkname[5]=='B') faktor=11; 
 if(brkname[5]=='C') faktor=12; 
 if(brkname[5]=='D') faktor=13; 
 if(brkname[5]=='E') faktor=14; 
 if(brkname[5]=='F') faktor=15; 
 if(brkname[5]=='G') faktor=16; 
 if(brkname[5]=='H') faktor=17; 
 if(brkname[5]=='I') faktor=18; 
 if(brkname[5]=='J') faktor=19; 
 if(brkname[5]=='K') faktor=20; 
 if(brkname[5]=='L') faktor=21; 
 if(brkname[5]=='M') faktor=22; 
 if(brkname[5]=='N') faktor=23; 
 if(brkname[5]=='O') faktor=24; 
 if(brkname[5]=='P') faktor=25; 
 if(brkname[5]=='Q') faktor=26; 
 if(brkname[5]=='R') faktor=27; 
 if(brkname[5]=='S') faktor=28; 
 if(brkname[5]=='T') faktor=29; 
 if(brkname[5]=='U') faktor=30; 
 if(brkname[5]=='V') faktor=31; 
 if(brkname[5]=='W') faktor=32; 
 if(brkname[5]=='X') faktor=33; 
 if(brkname[5]=='Y') faktor=34; 
 if(brkname[5]=='Z') faktor=35; 

 if(brkname[6]=='0') remaindr=0; 
 if(brkname[6]=='1') remaindr=1; 
 if(brkname[6]=='2') remaindr=2; 
 if(brkname[6]=='3') remaindr=3; 
 if(brkname[6]=='4') remaindr=4; 
 if(brkname[6]=='5') remaindr=5; 
 if(brkname[6]=='6') remaindr=6; 
 if(brkname[6]=='7') remaindr=7; 
 if(brkname[6]=='8') remaindr=8; 
 if(brkname[6]=='9') remaindr=9; 
 if(brkname[6]=='A') remaindr=10; 
 if(brkname[6]=='B') remaindr=11; 
 if(brkname[6]=='C') remaindr=12; 
 if(brkname[6]=='D') remaindr=13; 
 if(brkname[6]=='E') remaindr=14; 
 if(brkname[6]=='F') remaindr=15; 
 if(brkname[6]=='G') remaindr=16; 
 if(brkname[6]=='H') remaindr=17; 
 if(brkname[6]=='I') remaindr=18; 
 if(brkname[6]=='J') remaindr=19; 
 if(brkname[6]=='K') remaindr=20; 
 if(brkname[6]=='L') remaindr=21; 
 if(brkname[6]=='M') remaindr=22; 
 if(brkname[6]=='N') remaindr=23; 
 if(brkname[6]=='O') remaindr=24; 
 if(brkname[6]=='P') remaindr=25; 
 if(brkname[6]=='Q') remaindr=26; 
 if(brkname[6]=='R') remaindr=27; 
 if(brkname[6]=='S') remaindr=28; 
 if(brkname[6]=='T') remaindr=29; 
 if(brkname[6]=='U') remaindr=30; 
 if(brkname[6]=='V') remaindr=31; 
 if(brkname[6]=='W') remaindr=32; 
 if(brkname[6]=='X') remaindr=33; 
 if(brkname[6]=='Y') remaindr=34; 

 helper = faktor*35+remaindr;
 month = helper/100;
 minute = helper - month*100;

 if(month<10) {cmtname[5]='0'; sprintf(&cmtname[6],"%d",month);}
 if(month>9) {sprintf(monthnm,"%2d",month); strncpy(&cmtname[5],monthnm,2);}
 if(minute>9) {sprintf(minutnm,"%2d",minute); strncpy(&cmtname[11],minutnm,2);}
 if(minute<10) {cmtname[11]='0'; sprintf(&cmtname[12],"%d",minute);}
 
/* printf("%2d %2d %2d %2d\n",faktor,remaindr,month,minute); */

 if(brkname[3]=='1') {cmtname[7]='0';cmtname[8]='1';}
 if(brkname[3]=='2') {cmtname[7]='0';cmtname[8]='2';}
 if(brkname[3]=='3') {cmtname[7]='0';cmtname[8]='3';}
 if(brkname[3]=='4') {cmtname[7]='0';cmtname[8]='4';}
 if(brkname[3]=='5') {cmtname[7]='0';cmtname[8]='5';}
 if(brkname[3]=='6') {cmtname[7]='0';cmtname[8]='6';}
 if(brkname[3]=='7') {cmtname[7]='0';cmtname[8]='7';}
 if(brkname[3]=='8') {cmtname[7]='0';cmtname[8]='8';}
 if(brkname[3]=='9') {cmtname[7]='0';cmtname[8]='9';}
 if(brkname[3]=='A') {cmtname[7]='1';cmtname[8]='0';}
 if(brkname[3]=='B') {cmtname[7]='1';cmtname[8]='1';}
 if(brkname[3]=='C') {cmtname[7]='1';cmtname[8]='2';}
 if(brkname[3]=='D') {cmtname[7]='1';cmtname[8]='3';}
 if(brkname[3]=='E') {cmtname[7]='1';cmtname[8]='4';}
 if(brkname[3]=='F') {cmtname[7]='1';cmtname[8]='5';}
 if(brkname[3]=='G') {cmtname[7]='1';cmtname[8]='6';}
 if(brkname[3]=='H') {cmtname[7]='1';cmtname[8]='7';}
 if(brkname[3]=='I') {cmtname[7]='1';cmtname[8]='8';}
 if(brkname[3]=='J') {cmtname[7]='1';cmtname[8]='9';}
 if(brkname[3]=='K') {cmtname[7]='2';cmtname[8]='0';}
 if(brkname[3]=='L') {cmtname[7]='2';cmtname[8]='1';}
 if(brkname[3]=='M') {cmtname[7]='2';cmtname[8]='2';}
 if(brkname[3]=='N') {cmtname[7]='2';cmtname[8]='3';}
 if(brkname[3]=='O') {cmtname[7]='2';cmtname[8]='4';}
 if(brkname[3]=='P') {cmtname[7]='2';cmtname[8]='5';}
 if(brkname[3]=='Q') {cmtname[7]='2';cmtname[8]='6';}
 if(brkname[3]=='R') {cmtname[7]='2';cmtname[8]='7';}
 if(brkname[3]=='S') {cmtname[7]='2';cmtname[8]='8';}
 if(brkname[3]=='T') {cmtname[7]='2';cmtname[8]='9';}
 if(brkname[3]=='U') {cmtname[7]='3';cmtname[8]='0';}
 if(brkname[3]=='V') {cmtname[7]='3';cmtname[8]='1';}

 if(brkname[4]=='0') {cmtname[9]='0';cmtname[10]='0';}
 if(brkname[4]=='1') {cmtname[9]='0';cmtname[10]='1';}
 if(brkname[4]=='2') {cmtname[9]='0';cmtname[10]='2';}
 if(brkname[4]=='3') {cmtname[9]='0';cmtname[10]='3';}
 if(brkname[4]=='4') {cmtname[9]='0';cmtname[10]='4';}
 if(brkname[4]=='5') {cmtname[9]='0';cmtname[10]='5';}
 if(brkname[4]=='6') {cmtname[9]='0';cmtname[10]='6';}
 if(brkname[4]=='7') {cmtname[9]='0';cmtname[10]='7';}
 if(brkname[4]=='8') {cmtname[9]='0';cmtname[10]='8';}
 if(brkname[4]=='9') {cmtname[9]='0';cmtname[10]='9';}
 if(brkname[4]=='A') {cmtname[9]='1';cmtname[10]='0';}
 if(brkname[4]=='B') {cmtname[9]='1';cmtname[10]='1';}
 if(brkname[4]=='C') {cmtname[9]='1';cmtname[10]='2';}
 if(brkname[4]=='D') {cmtname[9]='1';cmtname[10]='3';}
 if(brkname[4]=='E') {cmtname[9]='1';cmtname[10]='4';}
 if(brkname[4]=='F') {cmtname[9]='1';cmtname[10]='5';}
 if(brkname[4]=='G') {cmtname[9]='1';cmtname[10]='6';}
 if(brkname[4]=='H') {cmtname[9]='1';cmtname[10]='7';}
 if(brkname[4]=='I') {cmtname[9]='1';cmtname[10]='8';}
 if(brkname[4]=='J') {cmtname[9]='1';cmtname[10]='9';}
 if(brkname[4]=='K') {cmtname[9]='2';cmtname[10]='0';}
 if(brkname[4]=='L') {cmtname[9]='2';cmtname[10]='1';}
 if(brkname[4]=='M') {cmtname[9]='2';cmtname[10]='2';}
 if(brkname[4]=='N') {cmtname[9]='2';cmtname[10]='3';}

 cmtname[13]=brkname[7];
 cmtname[14]='\0';

 printf("%s \n",cmtname);
 return 0; //K.Yuan

}
