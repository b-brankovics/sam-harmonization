FROM brankovics/snake-setup:v01

COPY src/ /opt/sam-harmonization

# Modify access rights so any user can run snakemake
RUN mkdir /.cache; chmod a+rwX /.cache

WORKDIR /analysis/snakemake

ENV PATH /opt/sam-harmonization:$PATH