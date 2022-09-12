{
  stdenv
, root
, roounfold
, cxxopts
, libyamlcpp
, boost
}:

stdenv.mkDerivation {
  pname = "misid-unfold";
  version = "0.2.6";

  src = builtins.path { path = ./../..; name = "misid-unfold"; };

  buildInputs = [ root roounfold cxxopts libyamlcpp boost ];

  installPhase = ''
    mkdir -p $out/bin
    cp bin/* $out/bin
  '';
}
