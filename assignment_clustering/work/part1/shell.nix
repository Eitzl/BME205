{ pkgs ? import <nixpkgs> {} }:
pkgs.mkShell {
  packages = [ (pkgs.python312.withPackages (p: [ p.numpy p.matplotlib ])) ];
}
