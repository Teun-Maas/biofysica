//
// mex_modbus.cc -- mex function that calls functions from libmodbus
// not all calls to libmodbus functions are implemented here (only the ones I needed until now),
// but adding calls is straightforward: define the subfunction and add it to the variable "lookupTable".
// 
// Written by Günter Windau
// 2014-04-22 version 1.0
//
#include <math.h>
#include <string.h>
#include <errno.h>
#include <matrix.h>
#include <mex.h>
#include <modbus/modbus.h>


static void mexModbus_new_tcp(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static void mexModbus_connect(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static void mexModbus_close(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static void mexModbus_free(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static void mexModbus_write_registers(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
static void mexModbus_read_registers(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);


typedef void (*pMexFunction)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

struct tLookupTable
{
   const char* name;
   const pMexFunction func;
};

tLookupTable lookupTable[] =
{
   { "new_tcp", mexModbus_new_tcp },
   { "connect", mexModbus_connect },
   { "close", mexModbus_close },
   { "free", mexModbus_free },
   { "write_registers", mexModbus_write_registers },
   { "read_registers", mexModbus_read_registers },
   { 0, 0 }
};

static const int MB_MAX_HANDLES = 5;
static modbus_t *ctx[MB_MAX_HANDLES] = { 0, };

void mexInit()
{
}

void mexExit()
{
   for (int i=0; i<MB_MAX_HANDLES; i++) {
      if (ctx[i]) {
         modbus_close(ctx[i]);
         modbus_free(ctx[i]);
      }
   }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   static int initialized = 0;

   if (!initialized) {
      mexAtExit(mexExit);
      mexInit();
      initialized = 1;
   }

   if (mxIsChar(prhs[0])) {
      char s[40];

      mxGetString(prhs[0], s, sizeof(s)-1);

      for (tLookupTable *p = lookupTable; p->name != 0; p++) {
         if (!strncasecmp(s, p->name, sizeof(s)-1)) {
            // name matched; call sub-function, discard the first rhs argument (i.e. the function name)
            p->func(nlhs, plhs, nrhs-1, prhs+1);
            return;
         }
      }
      mexErrMsgTxt("invalid or unimplemented function requested");
   }
   else {
      mexErrMsgTxt("no function requested");
   }
}

//
// argument checking functions
//
static int isScalar(const mxArray *p)
{
   int mrows = mxGetM(p);
   int ncols = mxGetN(p);
   if(!mxIsDouble(p) || mxIsComplex(p) || !(mrows==1 && ncols==1))
      return 0;
   return 1;
}

void convertArray(const mxArray *p, char **s)
{
  if (!mxIsChar(p))
      mexErrMsgTxt("argument must be a string.");

  int m = mxGetM(p);
  int n = mxGetN(p);
  mwSize buflen = m*n+1;
  *s = new char[buflen];
  if (mxGetString(p, *s, buflen)) {
     mexErrMsgTxt("string conversion error");
  }
}

void convertArray(const mxArray *p, int *i)
{
   int mrows = mxGetM(p);
   int ncols = mxGetN(p);

   if(!mxIsDouble(p) || mxIsComplex(p) || !(mrows==1 && ncols==1))
      mexErrMsgTxt("argument must be a scalar.");

   *i = mxGetScalar(p);
}

void convertArray(const mxArray *p, uint16_t **dest)
{
   int mrows = mxGetM(p);
   int ncols = mxGetN(p);

   if(!(mrows==1 || ncols==1))
      mexErrMsgTxt("argument must be a vector.");

   if (mxGetClassID(p) != mxUINT16_CLASS)
      mexErrMsgTxt("argument must be a vector of UINT16 values.");

   int n = mrows;         // choose either length
   if (n==1) n=ncols;

   *dest = (uint16_t *)mxGetData(p);
}

void checkNargs(int nargs, int min, int max, const char* errMsg)
{
   if (nargs<min || nargs>max)
      mexErrMsgTxt(errMsg);
}

void checkNrhs(int nargs, int min, int max)
{
   checkNargs(nargs, min, max, "incorrect number of input arguments");
}

void checkNlhs(int nargs, int min, int max)
{
   checkNargs(nargs, min, max, "incorrect number of output arguments");
}


static modbus_t *getCtx(int handle)
{
   if (handle >= MB_MAX_HANDLES) {
      mexErrMsgTxt("handle out of range");
   }
   if (ctx[handle]==0) {
      mexErrMsgTxt("invalid handle");
   }
   return ctx[handle];
}

static modbus_t *getCtx(const mxArray *p)
{
   int h;
   convertArray(p, &h);
   return getCtx(h);
}

static void fail(const char *msg = 0)
{
    const char *errMsg = modbus_strerror(errno);
    mexErrMsgTxt(errMsg);
}

// mexModbus_new_tcp(ip-address, port);
static void mexModbus_new_tcp(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   checkNlhs(nlhs, 0, 1);
   checkNrhs(nrhs, 2, 2);

   char *hostname;
   convertArray(prhs[0], &hostname);

   int port;
   convertArray(prhs[1], &port);

   for (int h=0; h<MB_MAX_HANDLES; h++) {
      if (ctx[h] == 0) {
         ctx[h] = modbus_new_tcp(hostname, port);
         delete hostname;
         if (ctx[h] == 0) {
            fail();
         }
         plhs[0]=mxCreateDoubleScalar(h); // return the new handle to the modbus context structure
         return;
      }
   }
   mexErrMsgTxt("out of modbus handles");
}

// modbus_connect(modbus_t *ctx)
static void mexModbus_connect(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   checkNlhs(nlhs, 0, 0);
   checkNrhs(nrhs, 1, 1);

   modbus_t *ctx = getCtx(prhs[0]);
   if (modbus_connect(ctx) < 0) {
      fail();
   }
}

// modbus_close(modbus_t *ctx);
static void mexModbus_close(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   checkNlhs(nlhs, 0, 0);
   checkNrhs(nrhs, 1, 1);

   modbus_t *ctx = getCtx(prhs[0]);
   if (ctx) {
      modbus_close(ctx);
   }
}

// modbus_free(modbus_t *ctx);
static void mexModbus_free(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   checkNlhs(nlhs, 0, 0);
   checkNrhs(nrhs, 1, 1);

   int h;
   convertArray(prhs[0], &h);
   if (ctx[h]) {
      modbus_free(ctx[h]);
      ctx[h]=0;
   }
}

// modbus_write_registers(modbus_t *ctx, int addr, int nb, const uint16_t *data);
static void mexModbus_write_registers(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   checkNlhs(nlhs, 0, 0);
   checkNrhs(nrhs, 3, 3);

   modbus_t *ctx = getCtx(prhs[0]);
   if (!ctx)
      fail();

   int addr;
   uint16_t *data;
   convertArray(prhs[1], &addr);
   convertArray(prhs[2], &data);
   int m = mxGetM(prhs[2]);
   int n = mxGetN(prhs[2]);
   int nb = m*n;

   if (modbus_write_registers(ctx, addr, nb, data) < 0) {
      fail();
   }
}

// modbus_read_registers(modbus_t *ctx, int addr, int nb, uint16_t *dest);
static void mexModbus_read_registers(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   checkNlhs(nlhs, 1, 1);
   checkNrhs(nrhs, 3, 3);

   modbus_t *ctx = getCtx(prhs[0]);

   int addr, nb;
   convertArray(prhs[1], &addr);
   convertArray(prhs[2], &nb);
   uint16_t data[nb];

   if (!ctx)
      fail();

   if (modbus_read_registers(ctx, addr, nb, data) < 0) {
      fail();
   }
   plhs[0] = mxCreateDoubleMatrix(nb, 1, mxREAL);
   double *r = mxGetPr(plhs[0]);
   for (int i=0; i<nb; i++) {
      r[i] = data[i];
   }
}


