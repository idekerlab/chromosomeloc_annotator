# docker/Makefile for chromosome_loc_annotator

# Name and tag for the Docker image
IMAGE_NAME ?= cytoscape-chromloc
IMAGE_TAG  ?= latest
IMAGE      := $(IMAGE_NAME):$(IMAGE_TAG)

# Optional registry, e.g. REGISTRY=myregistry.io/myorg
REGISTRY   ?=
PUSH_IMAGE := $(if $(REGISTRY),$(REGISTRY)/$(IMAGE_NAME):$(IMAGE_TAG),$(IMAGE))

# Build context is the app directory (parent of this docker dir)
BUILD_CONTEXT := ..

.PHONY: help build run shell push

help:
	@echo "Targets:"
	@echo "  build        Build the Docker image ($(IMAGE))"
	@echo "  run          Run image with default args (shows help)"
	@echo "  shell        Start an interactive shell in the container"
	@echo "  push         Push the image to a registry (set REGISTRY=...)"
	@echo ""
	@echo "Variables:"
	@echo "  IMAGE_NAME   (default: cytoscape-chromloc)"
	@echo "  IMAGE_TAG    (default: latest)"
	@echo "  REGISTRY     (default: empty; e.g. myregistry.io/myorg)"

build:
	docker build \
	  -f Dockerfile \
	  -t $(IMAGE) \
	  $(BUILD_CONTEXT)

run:
	docker run --rm $(IMAGE)

shell:
	docker run --rm -it \
	  --entrypoint /bin/bash \
	  $(IMAGE)

push:
ifeq ($(REGISTRY),)
	@echo "ERROR: REGISTRY is not set. Example:"
	@echo "  make push REGISTRY=myregistry.io/myorg"
	@exit 1
else
	docker tag $(IMAGE) $(PUSH_IMAGE)
	docker push $(PUSH_IMAGE)
	@echo "Pushed $(PUSH_IMAGE)"
endif
