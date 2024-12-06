import matplotlib.pyplot as plt
from matplotlib.sankey import Sankey

def draw_clean_energy_sankey():
    fig = plt.figure(figsize=(8.3, 6))
    ax = fig.add_subplot(1, 1, 1, xticks=[], yticks=[])
    
    # Define the flows and layout
    sankey = Sankey(ax=ax, scale=0.011, head_angle=90, format='%.0f', unit='%')

    # Add the flows with improved offsets and black outlines
    sankey.add(flows=[100, -1, -44, -9, -12, -34],  # Positive = forward flow, negative = losses
               labels=[None, None, None, None, None, None],  # Labels handled manually
               orientations=[0, 1, 1, 1, 1, 0],  # -1 = downward, 0 = rightward
               pathlengths=[0.5, 1., 1., 1., 1., 1.],  # Adjust lengths for separation
               facecolor="#ffcccc",  # Light pink
               edgecolor="black")   # Pure black outline

    # Render the Sankey diagram
    diagrams = sankey.finish()

    # Manually add text in boxes
    labels = [
        ("Chemical Power", "100%"),
        ("Incomplete\n Combustion", "1%"),
        ("Heat Loss\n (Heat)", "44%"),
        ("Heat Loss\n (Gas Power)", "9%"),
        ("KE Loss", "12%"),
        ("Thrust Power", "34%")
    ]

    y_pos = 0.7
    positions = [
        (-0.69, 0.),     # Chemical Power
        (0.4, y_pos),     # Incomplete Combustion
        (1., y_pos),   # Heat Loss (Heat)
        (1.4, y_pos),   # Heat Loss (Gas Power)
        (1.8, y_pos),   # Kinetic Energy Loss
        (2.6, -0.33)    # Thrust Power
    ]

    for (label, value), (x, y) in zip(labels, positions):
        ax.text(
            x, y, f"{label}\n{value}", 
            bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.3"),
            ha="center", va="center", fontsize=9
        )

    ax.title.set_fontsize(14)

    plt.tight_layout()
    plt.show()

# Run the function
draw_clean_energy_sankey()
