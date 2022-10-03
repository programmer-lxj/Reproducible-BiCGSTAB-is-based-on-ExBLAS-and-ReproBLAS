#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <ScalarVectors.h>
#include "ToolsMPI.h"

// #define PRINT_SEND_RESOTRNF_VECTORS 1

/*********************************************************************************/

void Synchonization (MPI_Comm Synch_Comm, char *message) {
	int my_id, i ; 

	MPI_Comm_rank(Synch_Comm, &my_id); 
	MPI_Barrier (Synch_Comm);
	printf ("(%d) %s\n", my_id, message);
	if (my_id == 0) printf ("Waiting ... \n");
	if (my_id == 0) scanf ("%d", &i);
	if (my_id == 0) printf (" ... done\n");
	MPI_Barrier (Synch_Comm);
}

/*********************************************************************************/

// Return true if the corresponding asynchonous communication,
// defined by data, has been finalized
int TestSimple (void *data) {
	int flag = 0;
	ptr_SimpleNode smpnode = (ptr_SimpleNode) data;
	
	// Verify if the communication has finalized
	MPI_Test (&(smpnode->req), &flag, &(smpnode->sta));
	if (flag) {
		// Remove the data included in the simple node
		MPI_Wait (&(smpnode->req), &(smpnode->sta));
		free (smpnode);
	}

	// Returns the result
	return flag;
}

// Return true if the corresponding asynchonous communication,
// defined by data, has been finalized
int TestPacket (void *data) {
	int flag = 0;
	ptr_PacketNode pcknode = (ptr_PacketNode) data;
	
	// Verify if the communication has finalized
	MPI_Test (&(pcknode->req), &flag, &(pcknode->sta));
	if (flag) {
		// Remove the data included in the pack
		MPI_Wait (&(pcknode->req), &(pcknode->sta));
		MPI_Type_free (&(pcknode->pack));
		free (pcknode);
	}

	// Returns the result
	return flag;
}

// Detect the lost messages whose destination is one process
// into the processes of communicator Err_Comm
void DetectErrorsMPI (MPI_Comm Err_Comm) {
	int my_id, flag= 0;
	MPI_Status st;

	// Definition of the variable my_id
	MPI_Comm_rank(Err_Comm, &my_id); 
	// Test if some message exists
	MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, Err_Comm, &flag, &st);
	if (flag) {
		printf ("%d --> (%d,%d)\n", my_id, st.MPI_SOURCE, st.MPI_TAG);
	}
}

/*********************************************************************************/

// Prepare the structures required to send/receive a SparseMatrix structure
// * spr refers to the SparseMatrix from where the data is obtained
// * size is the number of rows to be sent
// * weight is the number of nonzeros to be sent
// * pcknode, where the resulted packet appears
void MakeSprMatrixPacket (SparseMatrix spr, int size, int weight, ptr_PacketNode pcknode) {
	int k;
	int *lblq = pcknode->lblq;
	MPI_Aint *dspl = pcknode->dspl;
	MPI_Datatype *type = pcknode->type;
		
	// Definition of reference pointer
	pcknode->ptr = (unsigned char *) spr.vptr;
	// Definition of the required vectors to create the packet
	type[0] = MPI_INT   ; lblq[0] = size+1; dspl[0] = (MPI_Aint) spr.vptr;
	type[1] = MPI_INT   ; lblq[1] = weight; dspl[1] = (MPI_Aint) spr.vpos;
	type[2] = MPI_DOUBLE; lblq[2] = weight; dspl[2] = (MPI_Aint) spr.vval;
	type[3] = MPI_UB    ; lblq[3] = 1     ; dspl[3] = (MPI_Aint) (spr.vptr+size+1);
	for (k=3; k>=0; k--) dspl[k] -= dspl[0]; 
	// Creation of the packet
	MPI_Type_create_struct (4, lblq, dspl, type, &(pcknode->pack));
	MPI_Type_commit(&(pcknode->pack));
}

void MakeSprMatrixSendPacket (SparseMatrix spr, int *vlen, int dimL, int dspL, 
															ptr_PacketNode pcknode) {
	int k, weight, dspZ;
	int *lblq = pcknode->lblq;
	MPI_Aint *dspl = pcknode->dspl;
	MPI_Datatype *type = pcknode->type;
		
//	printf ("dimL = %d , dspL = %d\n", dimL, dspL);
//	PrintInts (vlen, spr.dim1);
//	PrintInts (spr.vptr, spr.dim1+1);
	// Definition of reference pointer
	pcknode->ptr = (unsigned char *) (vlen+dspL);
	// Definition of the required vectors to create the packet
	dspZ = spr.vptr[dspL]; weight = spr.vptr[dspL+dimL] - dspZ;
//	printf ("dspZ = %d , weight = %d\n", dspZ, weight);
	type[0] = MPI_INT   ; lblq[0] = dimL  ; dspl[0] = (MPI_Aint) (vlen+dspL     );
	type[1] = MPI_INT   ; lblq[1] = weight; dspl[1] = (MPI_Aint) (spr.vpos+dspZ );
	type[2] = MPI_DOUBLE; lblq[2] = weight; dspl[2] = (MPI_Aint) (spr.vval+dspZ );
	type[3] = MPI_UB    ; lblq[3] = 1     ; dspl[3] = (MPI_Aint) (vlen+dimL+dspL);
	for (k=3; k>=0; k--) dspl[k] -= dspl[0]; 
	// Creation of the packet
	MPI_Type_create_struct (4, lblq, dspl, type, &(pcknode->pack));
	MPI_Type_commit(&(pcknode->pack));
}

void MakeSprMatrixRecvPacket (SparseMatrix sprL, int nnzL, ptr_PacketNode pcknode) {
	int k, dimL = sprL.dim1;
	int *lblq = pcknode->lblq;
	MPI_Aint *dspl = pcknode->dspl;
	MPI_Datatype *type = pcknode->type;
		
//	printf ("nnzL = %d\n", nnzL);
	// Definition of reference pointer
	pcknode->ptr = (unsigned char *) (sprL.vptr+1);
	// Definition of the required vectors to create the packet
	type[0] = MPI_INT   ; lblq[0] = dimL; dspl[0] = (MPI_Aint) (sprL.vptr+1);
	type[1] = MPI_INT   ; lblq[1] = nnzL; dspl[1] = (MPI_Aint) sprL.vpos;
	type[2] = MPI_DOUBLE; lblq[2] = nnzL; dspl[2] = (MPI_Aint) sprL.vval;
	type[3] = MPI_UB    ; lblq[3] = 1   ; dspl[3] = (MPI_Aint) (sprL.vptr+1+dimL);
	for (k=3; k>=0; k--) dspl[k] -= dspl[0]; 
	// Creation of the packet
	MPI_Type_create_struct (4, lblq, dspl, type, &(pcknode->pack));
	MPI_Type_commit(&(pcknode->pack));
}

// Compute the number of nonzeros elements of a PermSprMatrixRecvPacket packet
// * prc_src is the processor from which the messages is sent
// * dimL is the number of rows to be received
// * comm is the communicator in which the messages is sent
int ComputeSprMatrixRecvWeights (int prc_src, int dimL, MPI_Comm comm) {
	int tam, tam_int, tam_double, tam_ub;
	MPI_Status sta;

	// Definition of sizes
	MPI_Type_size(MPI_INT, &tam_int);
	MPI_Type_size(MPI_DOUBLE, &tam_double);
	MPI_Type_size(MPI_UB, &tam_ub);
	MPI_Probe (prc_src, Tag_Send_Packet_Matrix_To_Leaf, comm, &sta);
	MPI_Get_count (&sta, MPI_BYTE, &tam);

	// Return the number of nonzeros included in a packet
	return (tam - (dimL*tam_int + tam_ub)) / (tam_int + tam_double);
}


int DistributeMatrix (SparseMatrix spr, int index, ptr_SparseMatrix sprL, int indexL,
												int *vdimL, int *vdspL, int root, MPI_Comm comm) {
	int myId, nProcs;
	int i, dim = spr.dim1, divL, rstL, dimL, dspL, nnzL;
	ptr_PacketNode pcknode;

	// Getiing the parameter of the communicator
	MPI_Comm_rank(comm, &myId); MPI_Comm_size(comm, &nProcs); 
	// Broadcasting the matrix dimension
	MPI_Bcast (&dim, 1, MPI_INT, root, MPI_COMM_WORLD); 

	// Calculating the vectors of sizes (vdimL) and displacements (vdspl)
	divL = (dim / nProcs); rstL = (dim % nProcs);
	for (i=0; i<nProcs; i++) vdimL[i] = divL + (i < rstL);
	vdspL[0] = 0; for (i=1; i<nProcs; i++) vdspL[i] = vdspL[i-1] + vdimL[i-1];
	dimL = vdimL[myId];	dspL = vdspL[myId];	
	
	// Distribution of the matrix, by blocks
	if (root == myId) {
		int *vlen = NULL;

		CreateInts (&vlen, dim); ComputeLengthfromHeader (spr.vptr, vlen, dim);
		for (i=0; i<nProcs; i++) {
			if (i != myId) {
				// Creating the message for each destination
				pcknode = (ptr_PacketNode) malloc (sizeof(PacketNode));
				MakeSprMatrixSendPacket (spr, vlen, vdimL[i], vdspL[i], pcknode);
				MPI_Send (pcknode->ptr, 1, pcknode->pack, i, Tag_Send_Packet_Matrix_To_Leaf, comm);
				MPI_Type_free (&(pcknode->pack));
				free (pcknode);
			}
		}
		nnzL = spr.vptr[dspL+dimL] - spr.vptr[dspL]; 
//		CreateSparseMatrix (sprL, indexL, dimL, dimL, nnzL, 0);
		CreateSparseMatrix (sprL, indexL, dimL, dim, nnzL, 0);
		CopyInts    (vlen+dspL, sprL->vptr+1, dimL);
		CopyInts    (spr.vpos+spr.vptr[dspL], sprL->vpos, nnzL);
		CopyDoubles (spr.vval+spr.vptr[dspL], sprL->vval, nnzL);

		RemoveInts (&vlen);
	} else {
		MPI_Status sta;

		// Compute the number of nonzeroes and creating the local matrix
		nnzL = ComputeSprMatrixRecvWeights (root, dimL, comm);
//		CreateSparseMatrix (sprL, indexL, dimL, dimL, nnzL, 0);
		CreateSparseMatrix (sprL, indexL, dimL, dim, nnzL, 0);
		// Receiving the data on the local matrix
		pcknode = (ptr_PacketNode) malloc (sizeof(PacketNode));
		MakeSprMatrixRecvPacket (*sprL, nnzL, pcknode);
		MPI_Recv (pcknode->ptr, 1, pcknode->pack, root, Tag_Send_Packet_Matrix_To_Leaf,
               comm, &sta);
		MPI_Type_free (&(pcknode->pack));
		free (pcknode);
	}
	*(sprL->vptr) = indexL; TransformLengthtoHeader (sprL->vptr, dimL);

	return dim;
}

