"""Render only observed local role events and executed MMFF conformers."""
import json
import numpy as np


def render(out, snapshot, events):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from swarm.memory_ledger import ROLES
    kinds = ("created", "claimed", "SUCCEEDED", "BLOCKED", "FAILED")
    counts = np.zeros((len(ROLES), len(kinds)), dtype=int)
    for event in events:
        body = json.loads(event["body"])
        role = body["payload"].get("role")
        if role in ROLES and body["kind"] in kinds:
            counts[ROLES.index(role), kinds.index(body["kind"])] += 1
    fig, ax = plt.subplots(figsize=(9, 4.8), layout="constrained")
    im = ax.imshow(counts, cmap="Greys", vmin=0, vmax=max(1, counts.max()), aspect="auto")
    ax.set_xticks(range(len(kinds)), kinds, rotation=20, ha="right", fontsize=9)
    ax.set_yticks(range(len(ROLES)), ROLES)
    for i in range(len(ROLES)):
        for j in range(len(kinds)):
            ax.text(j, i, str(counts[i,j]), ha="center", va="center",
                    color="white" if counts[i,j] > counts.max()/2 else "black")
    ax.set_title("Observed local role events (not LLM dialogue)\n"
                 f"Provider tokens: {snapshot['actual_provider_tokens']:,}; ceiling: {snapshot['token_cap']:,}")
    fig.colorbar(im, ax=ax, label="Recorded event count", ticks=sorted(set(counts.ravel())))
    fig.savefig(out / "figures/swarm_dialogue_dynamics.png", dpi=300)
    plt.close(fig)
    alpha = next(t for t in snapshot["tasks"] if t["role"] == "Alpha")
    result = json.loads(alpha["result"])
    if "artifact" not in result:
        return
    data = json.loads((out / result["artifact"]).read_text(encoding="utf-8"))
    fig, axes = plt.subplots(1, len(data["records"]), figsize=(11, 4.5), layout="constrained", sharey=True)
    all_y = [c["relative_energy_kcal_mol"] for r in data["records"] for c in r["conformers"] if c["converged"]]
    for ax, record in zip(np.atleast_1d(axes), data["records"]):
        entries = [e for e in record["conformers"] if e["converged"]]
        ax.scatter([e["ring_dihedral_deg"] for e in entries], [e["relative_energy_kcal_mol"] for e in entries],
                   marker="o", facecolors="none", edgecolors="black")
        ax.set_title(f"{record['id']}: {record['ring_size']}-member ring\n{len(entries)} converged conformers")
        ax.set_xlim(-180, 180)
        ax.set_ylim(-0.2, max(all_y, default=1) * 1.12 + 0.2)
        ax.set_xlabel("One ring dihedral (degrees)")
        ax.grid(alpha=.2)
    axes[0].set_ylabel("Within-molecule relative MMFF94s energy (kcal/mol)")
    fig.suptitle("Executed local lactam benchmarks — not absolute ring strain or binding free energy", fontsize=11)
    fig.savefig(out / "figures/macrocycle_strain_map.png", dpi=300)
    plt.close(fig)
