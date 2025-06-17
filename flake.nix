{
  description = "Unfolding Muon misID samples.";

  inputs = {
    root-curated.url = "github:umd-lhcb/root-curated/lmg_misid_validation_fit";
    nixpkgs.follows = "root-curated/nixpkgs";
    flake-utils.follows = "root-curated/flake-utils";
    pyTuplingUtils.url = "github:umd-lhcb/pyTuplingUtils";
    pidcalib2.url = "github:umd-lhcb/pidcalib2";
  };

  outputs = { self, nixpkgs, flake-utils, root-curated, pyTuplingUtils, pidcalib2 }:
    {
      overlay = import ./nix/overlay.nix;
    } //
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config = { allowUnfree = true; };
          overlays = [ root-curated.overlay pyTuplingUtils.overlay pidcalib2.overlay self.overlay ];
        };
        python = pkgs.python3;
        pythonPackages = python.pkgs;
      in
      rec {
        packages = flake-utils.lib.flattenTree {
          dev-shell = devShell.inputDerivation;
          inherit (pkgs)
            misid-unfold
            misid-unfold-applyer
          ;
        };
        devShell = pkgs.mkShell rec {
          name = "misid-unfold-dev";
          buildInputs = (with pkgs; with pythonPackages; [
            # Dev tools
            clang-tools

            root
            roounfold
            cxxopts
            libyamlcpp
            boost

            # Python stack
            pylint

            # Python dependencies
            statsmodels
            numpy
            pyyaml
            pythonPackages.pyTuplingUtils  # naming clash with the input!
            pythonPackages.pidcalib2

            # LaTeX
            (texlive.combine {
              inherit (texlive)
                scheme-basic
                # Explicit dependencies
                latexmk
                blkarray
                booktabs
                amsmath
                enumitem
                mathtools
                standalone
                preview
                # Implicit dependencies
                pgf
                pgfopts
                etoolbox
                translator
                selnolig
                everysel
                xkeyval
                ;
            })
          ]);

          FONTCONFIG_FILE = pkgs.makeFontsConf {
            fontDirectories = with pkgs; [
              gyre-fonts
            ];
          };

          shellHook = ''
            export PATH=$(pwd)/scripts:$(pwd)/bin:$PATH
            export MPLBACKEND=agg  # the backend w/o a UI
            export MPLCONFIGDIR=$(pwd)/.matplotlib
          '';
        };
      });
}
