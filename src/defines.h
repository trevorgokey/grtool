
#ifndef DEFINES_H
#define DEFINES_H

#define likely(x)   __builtin_expect(x, 1)
#define unlikely(x) __builtin_expect(x, 0)

#define true		  1
#define false		  0

#define TRAJECTORY_LOCATION   "wt.nc"

#define USE_DOUBLE          1
#define USE_RECYCLE         1
#define STATS 		        0
#define USE_CL 		        0


#define MAX_FRAMES 	        2
#define WORKER_QUEUE_SIZE   1
#define STRIDE 		        117
#define LAST 		        (50000)
#define CLEAN_UP            1
#define SINGLE_TX           1

//broken///
#define CL_ASYNC            0
#define USE_MAP 	        0
//////////

//useless?//
#define USE_PIN 	        0
///////////

#define READ_BLK            1024
#define TRANS_BLK           256
#define ALIGN_BLK           256

#define TRANS_SHR_SIZE      256
#define ALIGN_SHR_SIZE      256

#define TRANS               1
#define WRITE_TRANS         0
#define ALIGN               1
#define RMSD                0
#define WRITE               0

#define LIMIT_FRAME_CNT

#if USE_DOUBLE
typedef double data_size;
#else
typedef float data_size;
#endif

#endif /* DEFINES_H */
