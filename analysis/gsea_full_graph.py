import math
import matplotlib.pyplot as plt
from pathlib import Path


def create_seamless_aspect_grid(images, cols=6, save_path="grid_output.png"):
    total_images = len(images)
    rows = math.ceil(total_images / cols)
    
    # Target aspect ratio: width / height = 1882 / 1815
    aspect_ratio = 1882 / 1815  
    
    # Calculate figure dimensions so the grid layout fits naturally
    fig_w = 15  # Desired width in inches
    fig_h = fig_w * (rows / cols) / aspect_ratio
    
    # layout="compressed" packs subplots tightly when keeping aspect ratios fixed
    fig, axes = plt.subplots(rows, cols, figsize=(fig_w, fig_h), layout="compressed")
    
    axes = axes.flatten()
    
    for i in range(len(axes)):
        ax = axes[i]
        
        if i < total_images:
            img = plt.imread(images[i]) if isinstance(images[i], str) else images[i]
            
            ax.imshow(img, aspect='equal')
            ax.axis('off')  # Hide tick marks and frame
            
            # Ensure subplot box stays edge-to-edge
            ax.set_adjustable('box')
        else:
            # Delete remaining 2 empty slots (35 & 36) to prevent dead spacing
            fig.delaxes(ax)
            
    # Zero out any remaining subplot spacing margins
    fig.subplots_adjust(wspace=0, hspace=0, left=0, right=1, bottom=0, top=1)
    
    plt.savefig(save_path, dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show()

images_folder = Path("./results/gsea_results/prerank")

images_list = [str(img) for img in images_folder.glob("*.png")]

create_seamless_aspect_grid(images_list)