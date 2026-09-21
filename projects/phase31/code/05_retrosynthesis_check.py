"""Route-graph screening only, NOT a retrosynthesis service simulation."""
def audit_route(route):
    if route is None:
        return dict(status="UNKNOWN", synthesis_feasible=None, reason="No verified route supplied")
    nodes = route.get("steps", [])
    ids = [n["id"] for n in nodes]
    if not nodes or len(set(ids)) != len(ids):
        raise ValueError("unique nonempty reaction-step ids required")
    graph = {n["id"]: n.get("parents", []) for n in nodes}
    visiting, depth = set(), {}
    def visit(key):
        if key not in graph:
            raise ValueError("unknown precursor step")
        if key in visiting:
            raise ValueError("cyclic synthesis graph")
        if key in depth:
            return depth[key]
        visiting.add(key)
        depth[key] = 1 + max((visit(p) for p in graph[key]), default=0)
        visiting.remove(key)
        return depth[key]
    longest = max(visit(k) for k in graph)
    return dict(status="REJECT_POLICY" if longest > 12 else "STRUCTURAL_CHECK_ONLY",
                longest_linear_steps=longest, total_steps=len(nodes), synthesis_feasible=None,
                required=["reaction precedent", "stereochemistry", "ring closure selectivity", "yield and purification"])
