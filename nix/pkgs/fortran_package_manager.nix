{ lib, stdenv
, fetchurl
, fetchFromGitHub
, gfortran8
, git
, cacert
}:

stdenv.mkDerivation rec {
  pname = "fpm";
  version = "0.2.0";

  srcs = [
    (fetchurl {
    url = "https://github.com/fortran-lang/fpm/releases/download/v${version}/fpm-${version}.f90";
    sha256 = "sha256-edUEH1zr0a3/mZAXs7X4jYgUJn29D6qg98cgQR7eRj4=";
    })
    (fetchFromGitHub {
    owner = "fortran-lang";
    repo = "fpm";
    rev = "v${version}";
    sha256 = "sha256-bWIdDVNzw74hLYCzaSTMC6P5o+1nRqLQPYbtN9lN+Ps=";
  })
  ];

  nativeBuildInputs = [ gfortran8 git cacert ];

  unpackPhase = ''
  for _src in $srcs; do
    cp -r "$_src" $(stripHash "$_src")
  done
'';

  buildPhase = ''
  cp fpm-${version}.f90 source/fpm
  cd source/fpm
  gfortran -J . fpm-${version}.f90 -o fpm
'';

  installPhase = ''
  ./fpm install --flag "-g -fbacktrace -O3" --prefix $out
'';

  doCheck = false;

  meta = with lib; {
    description = "Fortran Package Manager (fpm)";
    homepage = "https://fpm.fortran-lang.org/";
    # platforms = platforms.unix ++ platforms.macos ++ platforms.windows;
    license = licenses.bsd3;
    maintainers = [ maintainers.HaoZeke ];
  };

}
