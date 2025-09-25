"""
EarthImager 2D Validation Data
Known reference values from toy-14-dd example for validation
"""

# Known values from your toy-14-dd_trial5 example
TOY_14DD_REFERENCE = {
    "survey": {
        "electrodes": 14,
        "measurements": 74,
        "electrode_spacing": 1.0,
        "array_type": "dipole-dipole"
    },
    "mesh_from_out_file": {
        "nodes": 468,
        "nodes_x": 39,
        "nodes_y": 12,
        "elements": 418,
        "elements_x": 38,
        "elements_y": 11,
        "parameters": 210,
        "parameters_x": 30,
        "parameters_y": 7
    },
    "expected_voltage_range": {
        "min": -7.42905,
        "max": -0.0985201
    },
    "expected_resistivity_range": {
        "min": 131.8,
        "max": 160.2
    },
    "three_layer_model": {
        "rho1": 100,  # ohm-m
        "h1": 0.45,   # m
        "rho2": 200,  # ohm-m  
        "h2": 1.03,   # m
        "rho3": 120   # ohm-m (basement)
    }
}

def validate_against_reference(processed_data: dict) -> dict:
    """Compare processed data against known reference values"""
    
    ref = TOY_14DD_REFERENCE
    validation_results = {}
    
    # Survey validation
    survey = processed_data.get("survey", {})
    validation_results["survey"] = {
        "electrodes": {
            "expected": ref["survey"]["electrodes"],
            "actual": survey.get("electrodes", 0),
            "match": survey.get("electrodes", 0) == ref["survey"]["electrodes"]
        },
        "measurements": {
            "expected": ref["survey"]["measurements"], 
            "actual": survey.get("measurements", 0),
            "match": survey.get("measurements", 0) == ref["survey"]["measurements"]
        },
        "spacing": {
            "expected": ref["survey"]["electrode_spacing"],
            "actual": survey.get("electrode_spacing", 0),
            "match": abs(survey.get("electrode_spacing", 0) - ref["survey"]["electrode_spacing"]) < 0.1
        }
    }
    
    # Mesh validation
    mesh = processed_data.get("mesh", {})
    validation_results["mesh"] = {
        "nodes_in_range": {
            "expected_range": "400-500",
            "actual": mesh.get("total_nodes", 0),
            "reasonable": 400 <= mesh.get("total_nodes", 0) <= 500
        },
        "nodes_x_in_range": {
            "expected_range": "35-45", 
            "actual": mesh.get("nodes_x", 0),
            "reasonable": 35 <= mesh.get("nodes_x", 0) <= 45
        },
        "parameters_in_range": {
            "expected_range": "200-250",
            "actual": mesh.get("parameters", 0),
            "reasonable": 200 <= mesh.get("parameters", 0) <= 250
        }
    }
    
    # Data range validation  
    data_ranges = processed_data.get("data_ranges", {})
    voltage_range = data_ranges.get("voltage_range", {})
    resistivity_range = data_ranges.get("resistivity_range", {})
    
    validation_results["data_ranges"] = {
        "voltage": {
            "expected_min": ref["expected_voltage_range"]["min"],
            "actual_min": voltage_range.get("min", 0),
            "min_reasonable": abs(voltage_range.get("min", 0) - ref["expected_voltage_range"]["min"]) < 2.0
        },
        "resistivity": {
            "expected_range": f"{ref['expected_resistivity_range']['min']}-{ref['expected_resistivity_range']['max']}",
            "actual_range": f"{resistivity_range.get('min', 0):.1f}-{resistivity_range.get('max', 0):.1f}",
            "range_reasonable": (
                120 <= resistivity_range.get("min", 0) <= 140 and
                150 <= resistivity_range.get("max", 0) <= 170
            )
        }
    }
    
    # Overall validation score
    all_checks = []
    for category in validation_results.values():
        for check in category.values():
            if isinstance(check, dict):
                all_checks.append(check.get("match", check.get("reasonable", False)))
    
    validation_results["overall"] = {
        "total_checks": len(all_checks),
        "passed_checks": sum(all_checks),
        "validation_score": sum(all_checks) / len(all_checks) if all_checks else 0,
        "status": "PASS" if sum(all_checks) / len(all_checks) > 0.8 else "FAIL"
    }
    
    return validation_results