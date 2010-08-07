C3E1v_anycolor.cg
!!ARBvp1.0
# cgc version 2.2.0006, build date Mar 31 2009
# command line args: -profile arbvp1
# source file: C3E1v_anycolor.cg
#vendor NVIDIA Corporation
#version 2.2.0.6
#profile arbvp1
#program main
#semantic main.MVP : state.matrix.mvp
#var float corner : $vin.ATTR2 : ATTR2 : 0 : 1
#var float4 data : $vin.ATTR3 : ATTR3 : 1 : 1
#var float2 data2 : $vin.ATTR4 : ATTR4 : 2 : 1
#var float4 data3 : $vin.ATTR5 : ATTR5 : 3 : 1
#var float4x4 MVP : state.matrix.mvp : c[1], 4 : 4 : 1
#var float4 main.position : $vout.POSITION : HPOS : -1 : 1
#var float4 main.color : $vout.COLOR : COL0 : -1 : 1
#var float2 main.texCoord : $vout.TEXCOORD0 : TEX0 : -1 : 1
#const c[0] = 0 1 -0.5 0.5
#const c[5] = 2 3 0.15915491 0.25
#const c[6] = 24.980801 -24.980801 -60.145809 60.145809
#const c[7] = 85.453789 -85.453789 -64.939346 64.939346
#const c[8] = 19.73921 -19.73921 -1 1
#const c[9] = -9 0.75
PARAM c[10] = { { 0, 1, -0.5, 0.5 },
		state.matrix.mvp,
		{ 2, 3, 0.15915491, 0.25 },
		{ 24.980801, -24.980801, -60.145809, 60.145809 },
		{ 85.453789, -85.453789, -64.939346, 64.939346 },
		{ 19.73921, -19.73921, -1, 1 },
		{ -9, 0.75 } };
TEMP R0;
TEMP R1;
TEMP R2;
TEMP R3;
TEMP R4;
TEMP R5;
ABS R0.x, vertex.attrib[2];
FLR R0.x, R0;
ADD R0.y, -R0.x, -R0.x;
SLT R0.z, vertex.attrib[2].x, c[0].x;
MAD R0.x, R0.y, R0.z, R0;
ADD R0.z, R0.x, -c[5].y;
ABS R1.y, R0.z;
ADD R0.y, R0.x, -c[0];
ADD R0.z, R0.x, -c[5].x;
ABS R0.y, R0;
ABS R0.z, R0;
SGE R0.x, c[0], R0.y;
SGE R1.x, c[0], R0.z;
MAD R0.xy, R0.x, c[0], c[0].z;
ADD R0.zw, -R0.xyxy, c[0].w;
MAD R0.xy, R0.zwzw, R1.x, R0;
ADD R0.zw, -R0.xyxy, c[0].xywz;
SGE R1.x, c[0], R1.y;
MAD R5.xy, R0.zwzw, R1.x, R0;
MUL R0.x, vertex.attrib[3].w, c[5].z;
FRC R1.w, R0.x;
MAD R0.y, vertex.attrib[3].w, c[5].z, -c[5].w;
FRC R0.w, R0.y;
ADD R1.xyz, -R0.w, c[0].xwyw;
MUL R2.xyz, R1, R1;
MAD R3.xyz, R2, c[6].xyxw, c[6].zwzw;
MAD R3.xyz, R3, R2, c[7].xyxw;
MAD R3.xyz, R3, R2, c[7].zwzw;
ADD R0.xyz, -R1.w, c[0].xwyw;
MUL R0.xyz, R0, R0;
MAD R1.xyz, R0, c[6].xyxw, c[6].zwzw;
MAD R1.xyz, R1, R0, c[7].xyxw;
MAD R1.xyz, R1, R0, c[7].zwzw;
MAD R1.xyz, R1, R0, c[8].xyxw;
MAD R3.xyz, R3, R2, c[8].xyxw;
MAD R1.xyz, R1, R0, c[8].zwzw;
SLT R4.x, R1.w, c[5].w;
SGE R4.yz, R1.w, c[9].xxyw;
MOV R0.xz, R4;
DP3 R0.y, R4, c[8].zwzw;
DP3 R0.xy, R1, -R0;
MAD R1.xyz, R3, R2, c[8].zwzw;
SLT R2.x, R0.w, c[5].w;
SGE R2.yz, R0.w, c[9].xxyw;
MOV R0.w, R0.y;
MOV R3.xz, R2;
DP3 R3.y, R2, c[8].zwzw;
DP3 R1.xy, R1, -R3;
MOV R0.z, -R1.x;
MOV R0.y, R1;
MUL R0.zw, R5.y, R0;
MAD R0.xy, R5.x, R0, R0.zwzw;
MUL R0.xy, R0, vertex.attrib[4];
ADD R1.xy, R5, c[0].w;
MOV R0.zw, c[0].xyxy;
ADD R0.xy, R0, vertex.attrib[3];
DP4 result.position.w, R0, c[4];
DP4 result.position.z, R0, c[3];
DP4 result.position.y, R0, c[2];
DP4 result.position.x, R0, c[1];
MAD result.texcoord[0].xy, R1, vertex.attrib[5].zwzw, vertex.attrib[5];
MOV result.color, c[0].y;
END
# 62 instructions, 6 R-regs
