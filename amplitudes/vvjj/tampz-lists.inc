c contains lists with numbers of amplitudes in madgraph output
c that are needed for comparison of hand-made code with mg result
 
	integer listavu(12),listzvu(12),listavl(12),listzvl(12),
     #		lista4u(16),listz4u(16)	
	

c lists of diagrams for comparison to madgraph:
c charged current (k=5,6):
c       vb contributions:
 
 	data listavu / 25,26,29,109,110,113,169,171,177,179,223,227 /
	data listzvu / 27,28,30,111,112,114,170,172,178,180,224,228 /
	
 	data listavl /67,68,71,151,152,155,173,175,181,183,225,229/
	data listzvl /69,70,72,153,154,156,174,176,182,184,226,230/
	
c	vll contributions:	
	data lista4u /17,19,21,23,35,37,39,41,
     #			101,103,105,107,119,121,123,125/
	data listz4u /18,20,22,24,36,38,40,42,
     #			102,104,108,106,120,122,124,126/
	
	
