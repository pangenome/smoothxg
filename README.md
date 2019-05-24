# xg0

## a handle-graphy reimplementation of the XG succinct graph index

This implements the core features of the XG index behind the PathHandleGraph API.
It is a partial fork of the XG index in vg, aiming to simplify the codebase and decouple it from protobuf based classes in vg.
In addition to this refactoring, the short term objective is to implement improvements that support read mapping performance by caching path positions in the node space of the graph, as such operations currently consume a large amount of runtime during short read mapping.
