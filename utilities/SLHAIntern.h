#define Slhadata(i) sd(2*(i-1)+im)
#define SlhaData(i) sd(2*(i-1)+1)

#define Decay(i) sd(2*OffsetDecays+i)
#define LengthDecay 2*LengthDecays

#define IdBits 10
#define Decay_Entry(id,n) ishft(id,IdBits)+n
#define Decay_Id(i) (iand(int(Decay(i)),-2**IdBits)/2**IdBits)
#define Decay_Next(i) ibits(int(Decay(i)),0,IdBits)

#include "SLHADefs.h"

