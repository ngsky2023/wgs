if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi
export PATH=$PATH:/share/work1/wangrr/local/simple/bin
# #gcc 6.3.0
export PATH=/share/public/software/gcc/gmp/bin:/share/public/software/gcc/mpfr/bin:/share/public/software/gcc/mpc/bin:/share/public/software/gcc/6.3.0/bin:$PATH
export LD_LIBRARY_PATH=/share/public/software/gcc/gmp/lib:/share/public/software/gcc/mpfr/lib:/share/public/software/gcc/mpc/lib:/share/public/software/gcc/6.3.0/lib64:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=/share/public/software/gcc/6.3.0/include:/share/public/software/gcc/gmp/include:/share/public/software/gcc/mpfr/include:/share/public/software/gcc/mpc/include:$C_INCLUDE_PATH

# zlib
export LD_LIBRARY_PATH=/share/public/software/zlib/1.2.11/lib:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=/share/public/software/zlib/1.2.11/include:$C_INCLUDE_PATH

# cmake 3.9.6
export PATH=/share/public/software/cmake/3.9.6/bin:$PATH

# ghc 8.2.2 cabal 2.2.0.0
export PATH=/share/public/software/ghc/8.2.2/bin:/share/public/software/cabal/bin:$PATH
export LD_LIBRARY_PATH=/share/public/software/ghc/8.2.2/lib:/share/public/software/cabal/lib/x86_64-linux-ghc-8.2.2:$LD_LIBRARY_PATH

# java 1.8
export JAVA_HOME=/share/public/software/java/1.8.0_162
export JRE_HOME=${JAVA_HOME}/jre
export CLASSPATH=.:${JAVA_HOME}/lib:${JRE_HOME}/lib
export PATH=${JAVA_HOME}/bin:$PATH

# R 3.4.0
export PATH=/share/public/software/pcre/8.40/bin:/share/public/software/xz/5.2.3/bin:/share/public/software/curl/7.59/bin:/share/public/software/R/3.4.0/bin:$PATH
export LD_LIBRARY_PATH=/share/public/software/pcre/8.40/lib:/share/public/software/xz/5.2.3/lib:/share/public/software/curl/7.59/lib:/share/public/software/R/3.4.0/lib:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=/share/public/software/pcre/8.40/include:/share/public/software/xz/5.2.3/include:/share/public/software/curl/7.59/include:$C_INCLUDE_PATH

# python 2.7.14
export PATH=/share/public/software/workspaces1/env1/bin/:$PATH
export LD_LIBRARY_PATH=/share/public/software/python/2.7.14/lib:$LD_LIBRARY_PATH

# python 3.6.3
export PATH=/share/public/software/python/3.6.3/bin:$PATH
export LD_LIBRARY_PATH=/share/public/software/python/3.6.3/lib:$LD_LIBRARY_PATH

# perl 5.24.1
export PATH=/share/public/software/perl/5.24.1/bin:$PATH
export PERL5LIB=/share/public/software/perl/5.24.1/lib/5.24.1:$PERL5LIB



export PATH=/share/public/software/samtools/1.4/:$PATH

PATH="/home/wangrr/perl5/bin${PATH:+:${PATH}}"; export PATH;
PERL5LIB="/home/wangrr/perl5/lib/perl5${PERL5LIB:+:${PERL5LIB}}"; export PERL5LIB;
PERL_LOCAL_LIB_ROOT="/home/wangrr/perl5${PERL_LOCAL_LIB_ROOT:+:${PERL_LOCAL_LIB_ROOT}}"; export PERL_LOCAL_LIB_ROOT;
PERL_MB_OPT="--install_base \"/home/wangrr/perl5\""; export PERL_MB_OPT;
PERL_MM_OPT="INSTALL_BASE=/home/wangrr/perl5"; export PERL_MM_OPT;



export PATH=/share/public/software/bwa/0.7.15/:$PATH
export PATH=/share/public/software/glpk/4.65/bin/:$PATH
export LD_LIBRARY_PATH=/share/public/software/glpk/4.65/lib:/share/public/software/htslib/1.4/lib/:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=/share/public/software/glpk/4.65/include:$C_INCLUDE_PATH

export PATH=/share/public/software/node/8.11.4/bin:$PATH

export PERL5LIB=/share/work1/wangrr/local/Plib/lib/perl5/
