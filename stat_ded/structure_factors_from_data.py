import mmtbx.maps.utils
from mmtbx import utils
import os
import sys

from libtbx.utils import Sorry

from iotbx import (
    reflection_file_reader,
    reflection_file_utils,
    crystal_symmetry_from_any,
    extract_xtal_data,
)
import iotbx.pdb

from cctbx import crystal

from scitbx.array_family import flex

# ------------------------------------------------------------------------------------------
def structure_factors_from_data(
    args, output_file_name, log=sys.stdout, use_output_directory=True
):
    """
    INPUTS  :  *args: files names

    OUTPUTS :
    """
    master_params = mmtbx.maps.maps_including_IO_master_params()
    processed_args = utils.process_command_line_args(
        args=args, log=log, master_params=master_params
    )
    working_phil = processed_args.params
    params = working_phil.extract()
    fmodel_data_file_format = params.maps.output.fmodel_data_file_format
    if (len(params.maps.map_coefficients) == 0) and (len(params.maps.map) == 0):
        working_phil = master_params.fetch(sources=[working_phil])
        params = working_phil.extract()
    if len(processed_args.pdb_file_names) > 1:
        raise Sorry("Only one model file is allowed as input.")
    if (params.maps.input.pdb_file_name is None) and (
        len(processed_args.pdb_file_names) == 1
    ):
        params.maps.input.pdb_file_name = processed_args.pdb_file_names[0]
    if not os.path.isfile(str(params.maps.input.pdb_file_name)):
        raise Sorry(
            "model file is not given: maps.input.pdb_file_name=%s is not a file"
            % str(params.maps.input.pdb_file_name)
        )
    if (
        (params.maps.input.reflection_data.file_name is None)
        and (params.maps.input.reflection_data.r_free_flags.file_name is None)
        and (len(processed_args.reflection_file_names) == 1)
    ):
        params.maps.input.reflection_data.file_name = (
            processed_args.reflection_file_names[0]
        )
    working_phil = master_params.format(python_object=params)
    pdb_inp = iotbx.pdb.input(file_name=params.maps.input.pdb_file_name)
    # get all crystal symmetries
    cs_from_coordinate_files = [pdb_inp.crystal_symmetry_from_cryst1()]
    cs_from_reflection_files = []

    for rfn in [
        params.maps.input.reflection_data.file_name,
        params.maps.input.reflection_data.r_free_flags.file_name,
    ]:
        if os.path.isfile(str(rfn)):
            try:
                cs_from_reflection_files.append(
                    crystal_symmetry_from_any.extract_from(rfn)
                )
            except KeyboardInterrupt:
                raise
            except RuntimeError:
                pass

    crystal_symmetry = None

    try:
        crystal_symmetry = crystal.select_crystal_symmetry(
            from_coordinate_files=cs_from_coordinate_files,
            from_reflection_files=cs_from_reflection_files,
        )
    except AssertionError as e:
        if "No unit cell and symmetry information supplied" in str(e):
            raise Sorry(
                "Missing or incomplete symmetry information.  This program "
                + "will only work with reflection file formats that contain both "
                + "unit cell and space group records, such as MTZ files."
            )
    reflection_files = []
    reflection_file_names = []
    for rfn in [
        params.maps.input.reflection_data.file_name,
        params.maps.input.reflection_data.r_free_flags.file_name,
    ]:
        if (os.path.isfile(str(rfn))) and (not rfn in reflection_file_names):
            reflection_files.append(
                reflection_file_reader.any_reflection_file(
                    file_name=rfn, ensure_read_access=False
                )
            )
            reflection_file_names.append(rfn)
    reflection_file_server = reflection_file_utils.reflection_file_server(
        crystal_symmetry=crystal_symmetry,
        force_symmetry=True,
        reflection_files=reflection_files,  # [],
        err=log,
    )
    #
    reflection_data_master_params = extract_xtal_data.data_and_flags_master_params(
        master_scope_name="reflection_data"
    )

    reflection_data_input_params = processed_args.params.get(
        "maps.input.reflection_data"
    )
    reflection_data_params = (
        reflection_data_master_params.fetch(reflection_data_input_params)
        .extract()
        .reflection_data
    )
    determine_data_and_flags_result = extract_xtal_data.run(
        reflection_file_server=reflection_file_server,
        parameters=reflection_data_params,
        keep_going=True,
    )

    f_obs = determine_data_and_flags_result.f_obs
    r_free_flags = determine_data_and_flags_result.r_free_flags
    # test_flag_value = determine_data_and_flags_result.test_flag_value
    if r_free_flags is None:
        r_free_flags = f_obs.array(data=flex.bool(f_obs.data().size(), False))
        # test_flag_value = None
    pdb_hierarchy = pdb_inp.construct_hierarchy(set_atom_i_seq=True)
    atom_selection_manager = pdb_hierarchy.atom_selection_cache()
    xray_structure = pdb_hierarchy.extract_xray_structure(
        crystal_symmetry=crystal_symmetry
    )
    # apply omit selection
    if params.maps.omit.selection is not None:
        omit_selection = atom_selection_manager.selection(
            string=params.maps.omit.selection
        )
        keep_selection = ~omit_selection
        xray_structure = xray_structure.select(selection=keep_selection)
        pdb_hierarchy = pdb_hierarchy.select(keep_selection)
        atom_selection_manager = pdb_hierarchy.atom_selection_cache()
    #

    utils.setup_scattering_dictionaries(
        scattering_table=params.maps.scattering_table,
        xray_structure=xray_structure,
        d_min=f_obs.d_min(),
        log=log,
    )
    if params.maps.wavelength is not None:
        if params.maps.scattering_table == "neutron":
            raise Sorry(
                "Wavelength parameter not supported when the neutron "
                + "scattering table is used."
            )
        xray_structure.set_inelastic_form_factors(
            photon=params.maps.wavelength, table="sasaki"
        )

    fmodel = utils.fmodel_simple(
        xray_structures=[xray_structure],
        scattering_table=params.maps.scattering_table,
        f_obs=f_obs,
        r_free_flags=r_free_flags,
        outliers_rejection=params.maps.input.reflection_data.outliers_rejection,
        skip_twin_detection=params.maps.skip_twin_detection,
        bulk_solvent_correction=params.maps.bulk_solvent_correction,
        anisotropic_scaling=params.maps.anisotropic_scaling,
    )
    # fmodel_info = fmodel.info()
    if (params.maps.output.directory is not None) and (use_output_directory):
        assert os.path.isdir(params.maps.output.directory)
        output_dir = params.maps.output.directory
    else:
        output_dir = os.getcwd()
    if params.maps.output.prefix is not None:
        file_name_base = os.path.join(
            output_dir, os.path.basename(params.maps.output.prefix)
        )
    else:
        file_name_base = params.maps.input.pdb_file_name
        if file_name_base.count(".") > 0:
            file_name_base = file_name_base[: file_name_base.index(".")]
    # if params.maps.output.include_r_free_flags:
    #     r_free_flags_output = fmodel.r_free_flags().average_bijvoet_mates()
    if params.maps.output.fmodel_data_file_format is not None:
        fmodel_file_name = (
            output_file_name + "_fmodel." + params.maps.output.fmodel_data_file_format
        )
        print(
            "Writing fmodel arrays (Fobs, Fcalc, m, ...) to %s file."
            % fmodel_file_name,
            file=log,
        )
        fmodel_file_object = open(fmodel_file_name, "w")
        fmodel.export(
            out=fmodel_file_object, format=params.maps.output.fmodel_data_file_format
        )
        fmodel_file_object.close()
    return fmodel_file_name