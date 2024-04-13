// Copyright (c) 2007 onward, Piotr Wendykier
// Copyright (c) 2024 Carl Reinke
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer. 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES；
// LOSS OF USE, DATA, OR PROFITS； OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <stdlib.h>
#include <fftw3.h>

long double randld()
{
  long double r = 0;
  long double d = RAND_MAX;
  while (1)
  {
    r += rand() / d;
    if (r + 1 / d == r)
      return r;
    d *= RAND_MAX;
  }
}

void generate(int n, fftwl_complex *in, fftwl_complex *out)
{
  for (int i = 0; i < n; ++i)
  {
    in[i][0] = -0.5L + randld();
    in[i][1] = -0.5L + randld();
  }

  fftwl_plan p = fftwl_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftwl_execute(p);
  fftwl_destroy_plan(p);  
}

void save(int n, const char *ext, fftwl_complex *cs)
{
  char s[100];
  snprintf(s, sizeof(s), "fftw.%d.%s", n, ext);

  FILE *f = fopen(s, "wb");
  for (int i = 0; i < n; ++i)
  {
    double c[2] = { cs[i][0], cs[i][1] };
    fwrite(c, sizeof(c), 1, f);
  }
  fclose(f);
}

int main(int argc, char **argv)
{
  srand(1415926535);

  while (!feof(stdin))
  {
    int n;
    fftwl_complex *in, *out;

    fscanf(stdin, "%i\n", &n);
    printf("n = %i\n", n);

    in = (fftwl_complex *)fftwl_malloc(sizeof(fftwl_complex) * n);
    out = (fftwl_complex *)fftwl_malloc(sizeof(fftwl_complex) * n);
    
    generate(n, in, out);

    save(n, "in", in);
    save(n, "out", out);

    fftwl_free(in);
    fftwl_free(out);
  }
}
