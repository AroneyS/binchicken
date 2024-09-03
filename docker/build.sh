
#!/bin/bash -eo pipefail

export BINCHICKEN_VERSION=$(binchicken --version)
export BINCHICKEN_DOCKER_VERSION=samuelaroney/binchicken:$BINCHICKEN_VERSION

cp ../binchicken.yml . && \
sed 's/BINCHICKEN_VERSION/'$BINCHICKEN_VERSION'/g' Dockerfile.in > Dockerfile && \
DOCKER_BUILDKIT=1 docker build -t $BINCHICKEN_DOCKER_VERSION . && \
  docker run $BINCHICKEN_DOCKER_VERSION coassemble --full-help | head && \
  echo "Seems good - now you just need to 'docker push $BINCHICKEN_DOCKER_VERSION'" && \
  echo "And then run AroneyS/binchicken-installation"
