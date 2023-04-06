version development

# When a WDL document is imported all of its structs are added to a global struct namespace.
# This enables structs to be used by their name alone, without the need for any "namespace." prefix.

struct RunEnv {
  String docker
  Int cpu
  Int memory
  Int disks
}
