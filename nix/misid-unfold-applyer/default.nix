{
  stdenv
, root
, cxxopts
, libyamlcpp
, boost
}:

stdenv.mkDerivation {
  pname = "misid-unfold-applyer";
  version = "0.2.5";

  src = builtins.path { path = ./../..; name = "misid-unfold"; };

  buildInputs = [ root cxxopts libyamlcpp boost ];

  buildPhase = "make applyer";

  installPhase = ''
    mkdir -p $out/bin
    cp bin/* $out/bin
  '';
}
