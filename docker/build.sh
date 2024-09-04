
#!/bin/bash -eo pipefail

export BINCHICKEN_VERSION=$(binchicken --version)
export BINCHICKEN_DOCKER_VERSION=samuelaroney/binchicken:$BINCHICKEN_VERSION

cp ../binchicken.yml . && \
sed 's/BINCHICKEN_VERSION/'$BINCHICKEN_VERSION'/g' Dockerfile.in > Dockerfile && \
DOCKER_BUILDKIT=1 docker build -t $BINCHICKEN_DOCKER_VERSION . && \
  docker run $BINCHICKEN_DOCKER_VERSION coassemble --full-help | head && \
  docker run -v `pwd`/..:/data $BINCHICKEN_DOCKER_VERSION coassemble --forward test/data/sample_1.1.fq test/data/sample_2.1.fq test/data/sample_3.1.fq --reverse test/data/sample_1.2.fq test/data/sample_2.2.fq test/data/sample_3.2.fq --genomes test/data/GB_GCA_013286235.1.fna --singlem-metapackage test/data/singlem_metapackage.smpkg --assemble-unmapped --unmapping-max-identity 99 --unmapping-max-alignment 90 --prodigal-meta --output /data/test_docker && \
  echo "Seems good - now you just need to 'docker push $BINCHICKEN_DOCKER_VERSION'" && \
  echo "And then run AroneyS/binchicken-installation"
