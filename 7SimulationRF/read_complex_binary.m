function v = read_complex_binary (filename, count, startind)

  %% usage: read_complex_binary (filename, [count])
  %%
  %%  open filename and return the contents as a column vector,
  %%  treating them as 32 bit complex numbers
  %%

  floatSize = 8;
  m = nargchk (1,3,nargin);
  if (m)
    usage (m);
  end

  if (nargin < 2)
    count = Inf;
  end
  if (nargin < 3)
      startind = 0;
  end

  f = fopen (filename, 'rb');
  if (f < 0)
    v = 0;
  else
    seekWorked = fseek(f,floatSize*startind,-1);
    if seekWorked==-1
        fclose (f);
        error('Setting start index for read_complex_binary did not work')
        return
    end
    t = fread (f, [2, count], 'float');
    fclose (f);
    v = t(1,:) + t(2,:)*i;
    [r, c] = size (v);
    v = reshape (v, c, r);
  end
