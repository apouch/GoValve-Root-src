import os

class GreedyHelper:
    def __init__(self, location):
        self.greedy = location

    def run_reg(self, img_fix, img_mov, regout_deform_inv, mask_fix, \
        affine_init = '', regout_affine = '', regout_deform = '', reference_image = '', \
        multi_res_schedule = '100x30', metric_spec = 'SSD', threads = -1):

        """
        Make system call to greedy

        Calls greedy for deformable registration between two images. Creates
        affine transformation file and/or deformation field depending on which
        optional input filenames are provided.

        INPUT:
        - img_fix: fixed (reference) image filename
        - img_mov: moving image filename
        - affine_init: filename of affine transform for initialization (optional)
        - regout_affine: filename of output affine transformation (optional)
        - regout_deform: filename of output deformation (optional)
        - regout_deform_inv: filename of output inverse deformation (optional)
        - mask_fix: filename of mask for the fixed image
        
        The optional input arguments determine which parameters are used in the 
        call to greedy (affine, deformable registration, or both). 

        Use full paths for all filenames.
        """

        cmdbase = f'{self.greedy} -d 3 '
        if threads > 0:
            cmdbase = cmdbase + f'-threads {threads} '

        if regout_affine != '' and regout_deform != '':
            # Affine generation
            aff_cmd = cmdbase + f'-a -i {img_fix} {img_mov} '

            if reference_image != '':
                aff_cmd = f'{aff_cmd} -rf {reference_image} '
            
            aff_cmd = f'{aff_cmd} \
                -ia-identity \
                -dof 6 \
                -s 3mm 1.5mm \
                -gm {mask_fix} \
                -o {regout_affine} '
            
            print('greedy_call (affine): ', aff_cmd)
            os.system(aff_cmd)

            # Deform generation
            def_cmd = f'{cmdbase} -i {img_fix} {img_mov} -it {regout_affine} '

            if reference_image != '':
                def_cmd = f'{def_cmd} -rf {reference_image} '

            def_cmd = f'{def_cmd} \
                -m {metric_spec} \
                -n {multi_res_schedule} \
                -s 3mm 1.5mm \
                -gm {mask_fix} \
                -o {regout_deform} '

            if regout_deform_inv != '':
                def_cmd = f'{def_cmd} -oinv {regout_deform_inv} '

            print('greedy_call (deformable): ', def_cmd)
            os.system(def_cmd)
            #print('greedy_call: Affine + Deformable registration computed!')

        elif regout_affine == '' and regout_deform != '':
            if regout_deform_inv != '':
                if affine_init != '':
                    cmd = f'{cmdbase} -i {img_fix} {img_mov} '

                    if reference_image != '':
                        cmd = f'{cmd} -rf {reference_image} '
                    
                    cmd = f'{cmd} -m {metric_spec} \
                        -n {multi_res_schedule} \
                        -it {affine_init} \
                        -gm {mask_fix} \
                        -s 3mm 1.5mm \
                        -o {regout_deform} \
                        -oinv {regout_deform_inv}'
                    #print('Initialized with input affine transform')
                    print('greedy_call: ', cmd)
                    os.system(cmd)
                    #print('greedy_call: Only deformable registration computed!')
        
                
    def apply_warp(self, image_type, img_fix, img_mov, img_reslice, \
        reg_affine = '', reg_deform = '', reg_deform_inv = '', threads = -1):
        """
        Calls greedy to apply an affine and/or other deformation to a grayscale
        image, label map (segmentation), or vtk mesh. 

        INPUT:
        - image_type: 'grayscale', 'label', or 'mesh'
        - img_fix: fixed (reference) image filename
        - img_reslice: filename of output (warped) image
        - reg_affine: filename of affine registration (optional)
        - reg_deform: filename of deformation field (optional)
        - reg_deform_inv: filename of inverse deformation field (optional)
        
        The optional input arguments determine which parameters are used in the 
        call to greedy (affine, deformable registration, or both). 
        
        Use full paths for all filenames.
        """
        cmdbase = f'{self.greedy} -d 3 '
        if threads > 0:
            cmdbase = cmdbase + f'-threads {threads} '

        cmd = f'{cmdbase} -d 3 \
            -rf {img_fix} \
            -r {reg_deform} {reg_affine} '
        
        if image_type == 'grayscale':
            cmd = cmd + f' \
                -ri LINEAR \
                -rm {img_mov} {img_reslice} '
        elif image_type == 'label':
            cmd = cmd + f' \
                -ri NN \
                -rm {img_mov} {img_reslice} '
        elif image_type == 'mesh':
            cmd = cmd + f' \
                -rs {img_mov} {img_reslice} '

        print('apply_warp: ', cmd)
        os.system(cmd)
        #print('apply_warp: Affine + Deformable transformation applied!')
