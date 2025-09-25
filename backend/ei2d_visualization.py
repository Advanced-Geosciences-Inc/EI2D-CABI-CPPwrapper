"""
EarthImager 2D Visualization Service
Integrates the OUT file plotting capabilities with the web interface
"""

import os
import sys
import tempfile
import base64
import io
from pathlib import Path
from typing import Dict, List, Any, Optional
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for web
import matplotlib.pyplot as plt
import numpy as np

# Add the EarthImager directory to Python path to import the plotting code
sys.path.append('/app/earthimager')

class EI2DVisualizationService:
    """Web service wrapper for EarthImager 2D plotting functionality"""
    
    def __init__(self):
        self.temp_dir = Path(tempfile.mkdtemp(prefix="ei2d_plots_"))
        
    def generate_plots_from_out_content(self, out_file_content: str, 
                                      plot_options: Optional[Dict] = None) -> Dict[str, Any]:
        """Generate all EarthImager 2D plots from OUT file content"""
        
        try:
            # Create temporary OUT file
            out_file_path = self.temp_dir / "results.out"
            out_file_path.write_text(out_file_content)
            
            # Import the plotting class
            from generate_all_plots import ERTDataProcessor
            
            # Set default plot options
            options = {
                "inversion_slope_factor": 0.5,
                "pseudosection_slope_factor": 0.3,
                "all_iterations_pseudosections": False,
                "colormap": "jet",
                "show_contours": False
            }
            if plot_options:
                options.update(plot_options)
            
            # Create processor with temp output directory
            plot_output_dir = self.temp_dir / "plots"
            processor = ERTDataProcessor(
                out_file_path=str(out_file_path),
                output_dir=str(plot_output_dir),
                target_max_depth=options.get("target_max_depth"),
                bottom_width=options.get("bottom_width", 0.0)
            )
            
            # Generate all plots
            processor.generate_all_plots(
                inversion_slope_factor=options["inversion_slope_factor"],
                pseudosection_slope_factor=options["pseudosection_slope_factor"],
                all_iterations_pseudosections=options["all_iterations_pseudosections"],
                colormap=options["colormap"],
                show_contours=options["show_contours"]
            )
            
            # Collect generated plots as base64 images
            plots = self._collect_generated_plots(plot_output_dir)
            
            # Extract summary statistics
            summary = self._extract_plot_summary(processor)
            
            return {
                "success": True,
                "plots": plots,
                "summary": summary,
                "plot_options": options,
                "message": f"Generated {len(plots)} plots successfully"
            }
            
        except Exception as e:
            return {
                "success": False,
                "error": str(e),
                "message": "Plot generation failed"
            }
    
    def _collect_generated_plots(self, plot_dir: Path) -> Dict[str, str]:
        """Collect all generated plot images as base64 strings"""
        plots = {}
        
        if not plot_dir.exists():
            return plots
        
        # Define expected plot categories and their file patterns
        plot_categories = {
            "inversion_results": ["*inverted_resistivity*.png", "*final_iteration*.png"],
            "measured_data": ["*measured_apparent_resistivity*.png"],
            "calculated_data": ["*calculated_apparent_resistivity*.png"],
            "convergence": ["*convergence*.png"],
            "data_fit": ["*crossplot*.png", "*misfit*.png"],
            "all_iterations": ["*all_iterations*.png"]
        }
        
        for category, patterns in plot_categories.items():
            category_plots = []
            for pattern in patterns:
                for plot_file in plot_dir.glob(pattern):
                    try:
                        # Convert image to base64
                        with open(plot_file, 'rb') as f:
                            img_data = base64.b64encode(f.read()).decode('utf-8')
                            category_plots.append({
                                "filename": plot_file.name,
                                "image_data": img_data,
                                "title": self._generate_plot_title(plot_file.name)
                            })
                    except Exception as e:
                        print(f"Error processing {plot_file}: {e}")
            
            if category_plots:
                plots[category] = category_plots
        
        return plots
    
    def _generate_plot_title(self, filename: str) -> str:
        """Generate human-readable title from filename"""
        title = filename.replace(".png", "").replace("_", " ").title()
        
        # Clean up common abbreviations
        replacements = {
            "App Res": "Apparent Resistivity",
            "Vi": "V/I",
            "Calc": "Calculated", 
            "Meas": "Measured",
            "Iter": "Iteration"
        }
        
        for old, new in replacements.items():
            title = title.replace(old, new)
            
        return title
    
    def _extract_plot_summary(self, processor) -> Dict[str, Any]:
        """Extract summary information from the processor"""
        
        summary = {
            "mesh_info": processor.mesh_info,
            "electrode_spacing": processor.electrode_spacing,
            "num_electrodes": len(processor.electrode_positions),
            "num_measurements": len(processor.measured_data),
            "convergence_iterations": len(processor.convergence_data),
            "array_types": []
        }
        
        # Extract array types
        if processor.measured_data:
            array_types = set(d.get('array_type', 'unknown') for d in processor.measured_data)
            summary["array_types"] = list(array_types)
        
        # Extract resistivity ranges
        if processor.resistivity_models_per_iter:
            final_model = processor.resistivity_models_per_iter[-1]
            summary["final_resistivity_range"] = {
                "min": float(np.min(final_model)),
                "max": float(np.max(final_model)),
                "mean": float(np.mean(final_model)),
                "std": float(np.std(final_model))
            }
        
        # Extract convergence info
        if processor.convergence_data:
            final_convergence = processor.convergence_data[-1]
            summary["convergence_info"] = {
                "final_iteration": final_convergence.get("iteration", 0),
                "final_rms": final_convergence.get("rms_error", 0),
                "total_iterations": len(processor.convergence_data)
            }
        
        return summary
    
    def generate_specific_plot(self, out_file_content: str, plot_type: str, 
                              options: Optional[Dict] = None) -> Dict[str, Any]:
        """Generate a specific type of plot"""
        
        try:
            # Create temporary OUT file
            out_file_path = self.temp_dir / "results.out"
            out_file_path.write_text(out_file_content)
            
            from generate_all_plots import ERTDataProcessor
            
            plot_output_dir = self.temp_dir / "single_plot"
            processor = ERTDataProcessor(
                out_file_path=str(out_file_path),
                output_dir=str(plot_output_dir)
            )
            
            # Generate specific plot based on type
            plt.figure(figsize=(12, 8))
            
            if plot_type == "convergence":
                processor.plot_convergence_curve()
            elif plot_type == "measured_pseudosection":
                processor.plot_measured_apparent_resistivity()
            elif plot_type == "calculated_pseudosection":
                processor.plot_calculated_apparent_resistivity()
            elif plot_type == "inversion_result":
                processor.plot_and_export_inverted_sections_all_iterations()
            elif plot_type == "crossplot":
                processor.plot_crossplot()
            else:
                raise ValueError(f"Unknown plot type: {plot_type}")
            
            # Save to buffer
            buf = io.BytesIO()
            plt.savefig(buf, format='png', dpi=150, bbox_inches='tight')
            buf.seek(0)
            
            # Convert to base64
            img_data = base64.b64encode(buf.read()).decode('utf-8')
            plt.close()
            
            return {
                "success": True,
                "plot_data": img_data,
                "plot_type": plot_type,
                "message": f"Generated {plot_type} plot successfully"
            }
            
        except Exception as e:
            return {
                "success": False,
                "error": str(e),
                "message": f"Failed to generate {plot_type} plot"
            }
    
    def cleanup(self):
        """Clean up temporary files"""
        try:
            import shutil
            shutil.rmtree(self.temp_dir)
        except:
            pass


class PlotOptionsValidator:
    """Validate and provide defaults for plot options"""
    
    @staticmethod
    def validate_options(options: Dict) -> Dict:
        """Validate and set default values for plot options"""
        
        validated = {
            "inversion_slope_factor": float(options.get("inversion_slope_factor", 0.5)),
            "pseudosection_slope_factor": float(options.get("pseudosection_slope_factor", 0.3)),
            "all_iterations_pseudosections": bool(options.get("all_iterations_pseudosections", False)),
            "colormap": str(options.get("colormap", "jet")),
            "show_contours": bool(options.get("show_contours", False)),
            "target_max_depth": options.get("target_max_depth"),
            "bottom_width": float(options.get("bottom_width", 0.0))
        }
        
        # Validate colormap
        valid_colormaps = ['cividis', 'viridis', 'inferno', 'coolwarm', 'grayscale', 'jet']
        if validated["colormap"] not in valid_colormaps:
            validated["colormap"] = "jet"
        
        # Validate slope factors (0.1 to 1.0)
        validated["inversion_slope_factor"] = max(0.1, min(1.0, validated["inversion_slope_factor"]))
        validated["pseudosection_slope_factor"] = max(0.1, min(1.0, validated["pseudosection_slope_factor"]))
        
        return validated


# Global instance for the service
visualization_service = EI2DVisualizationService()