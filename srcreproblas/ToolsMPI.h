#ifndef ToolsMPI

#define ToolsMPI 1

// #include <SparseMatricesNew.h>
#include <SparseProduct.h>

/*********************************************************************************/

#define Tag_Demand_Matrix_From_Root       1001
#define Tag_Send_Task_To_Leaf             1002
#define Tag_Receive_Dims_Factor_From_Leaf 1003
#define Tag_End_Distribution_To_Leaf      1004
#define Tag_Send_Dims_Matrix_To_Leaf      1006
#define Tag_Send_Data_Matrix_To_Leaf      1007

#define Tag_Demand_Vector_From_Root       1011
#define Tag_Send_Dims_Vector_To_Father    1015
#define Tag_Send_Data_Vector_To_Father    1016

#define Tag_Send_Task_To_Root             1021
#define Tag_Send_Solution_To_Root         1022
#define Tag_Send_Dims_Vector_To_Children  1025
#define Tag_Send_Data_Vector_To_Children  1026

#define Tag_End_Resolution_To_Leaf        1031

#define Tag_Send_Vector_Up_1              1041
#define Tag_Send_Vector_Up_2              1042

#define Tag_Send_Packet_Matrix_To_Leaf     210
#define Tag_Receive_Data_Factor_From_Leaf  220
#define Tag_Send_Vector_To_Leaf            230
// #define Tag_FactorVector           240

/*********************************************************************************/

// typedef struct SimpleNode {
typedef struct {
	MPI_Status sta;
	MPI_Request req;
} SimpleNode, *ptr_SimpleNode;

// Return true if the corresponding asynchonous communication,
// defined by data, has been finalized
extern int TestSimple (void *data);

/*********************************************************************************/

// #define MaxPacketSize                    10000
#define MaxPacketSize                    5000

// typedef struct PacketNode {
typedef struct {
	unsigned char *ptr;
	int lblq[2*MaxPacketSize+3], vlen[MaxPacketSize];
	MPI_Aint dspl[2*MaxPacketSize+3];
	MPI_Datatype type[2*MaxPacketSize+3];
	MPI_Datatype pack;
	MPI_Status sta;
	MPI_Request req;
} PacketNode, *ptr_PacketNode;


/*********************************************************************************/

extern void Synchonization (MPI_Comm Synch_Comm, char *message);

/*********************************************************************************/

// Return true if the corresponding asynchonous communication,
// defined by data, has been finalized
extern int TestSimple (void *data);

// Return true if the corresponding asynchonous communication,
// defined by data, has been finalized
extern int TestPacket (void *data);

// Detect the lost messages whose destination is one process
// into the processes of communicator Err_Comm
extern void DetectErrorsMPI (MPI_Comm Err_Comm);

/*********************************************************************************/

// Prepare the structures required to send/receive a SparseMatrix structure
// * spr refers to the SparseMatrix from where the data is obtained
// * size is the number of rows to be sent
// * weight is the number of nonzeros to be sent
// * pcknode, where the resulted packet appears
extern void MakeSprMatrixPacket (SparseMatrix spr, int size, int weight, ptr_PacketNode pcknode);

extern void MakeSprMatrixSendPacket (SparseMatrix spr, int *len, int dimL, int dspL, 
																			ptr_PacketNode pcknode);

extern void MakeSprMatrixRecvPacket (SparseMatrix sprL, int nnzL, ptr_PacketNode pcknode);

// Compute the number of nonzeros elements of a PermSprMatrixRecvPacket packet
// * prc_src is the processor from which the messages is sent
// * sizes is the number of rows to be received
// * comm is the communicator in which the messages is sent
extern int ComputeSprMatrixRecvWeights (int prc_src, int sizes, MPI_Comm comm);

extern int DistributeMatrix (SparseMatrix spr, int index, ptr_SparseMatrix sprL, int indexL,
															int *vdimL, int *vdspL, int root, MPI_Comm comm);

#endif
