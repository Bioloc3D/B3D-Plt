import customtkinter as ctk
from CTkMessagebox import CTkMessagebox
from PIL import Image
from customtkinter import filedialog
import subprocess
import os
import platform
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import dabest
import warnings
import sys
import threading

matplotlib.use('pdf')
warnings.simplefilter("ignore")

def resource_path(relative_path):
    """Retourne le chemin absolu du fichier, compatible avec PyInstaller"""
    try:
        base_path = sys._MEIPASS  # Dossier temporaire utilisé par PyInstaller
    except AttributeError:
        base_path = os.path.abspath(".")

    full_path = os.path.join(base_path, relative_path)
    return full_path

logo_png = resource_path('/ressource/logo_png.png')

def quantification_script(pathway_to_raw_data, nb_cond, cond_names, save_path, colors=None):
    condition = {}
    cond = {}

    if not pathway_to_raw_data.endswith('/'):
        pathway_to_raw_data += '/'

    for i in ['C1 plots', 'C2 plots', 'C1xCo plots', 'C2xCo plots', 'Coloc plots']:
        if not os.path.isdir(f'{save_path}/{i}'):
            os.makedirs(f'{save_path}/{i}')

    C1_save_path = f'{save_path}/C1 plots'
    C2_save_path = f'{save_path}/C2 plots'
    Coloc_save_path = f'{save_path}/Coloc plots'
    C1xCo_save_path = f'{save_path}/C1xCo plots'
    C2xCo_save_path = f'{save_path}/C2xCo plots'

    for i in range(1, nb_cond + 1):
        condition[i] = {'nb': str(i), 'cond_name': cond_names[i - 1]}
        cond[i] = {}
    box_pairs = []
    sns.set_theme(style="ticks", palette="deep")

    # List preparation
    data_to_plot = {}
    for i in range(1, 39):
        data_to_plot[i] = pd.Series(dtype='object')

    name_list = []

    index = 0

    for cond in range(1, nb_cond + 1):
        for name in sorted(glob.glob(f'{pathway_to_raw_data}Cond{cond}_*')):
            # print(name) # Debugging
            index += 1
            # Open Excel file
            data = pd.read_excel(name,
                                 sheet_name=None,  # Takes all spreadsheet and store as a dict
                                 header=1,  # Second line is used as header
                                 )
            if 'C1' and 'C2' in data:  # Only keep dataset with C1 and C2 recorded elements
                name_list.append(name)
                # Find volume of stack
                volume_stack = data['Stack_']['Volume']

                # Find nucleus information if analysed
                if 'Nb of soma' in data['C1'].columns:
                    nucleus_number_C1 = data['C1']['Nb of soma'][0]
                    nucleus_mean_area_C1 = data['C1']['Mean area micron^2'][0]

                    nucleus_number_C2 = data['C2']['Nb of soma'][0]
                    nucleus_mean_area_C2 = data['C2']['Mean area micron^2'][0]

                # Find mean volume, total volume and number of C1 elements
                C1_elements = data['C1']['Total C1'][0]
                if data['C1'].shape[0] > 1:
                    if data['C1']['Label'][0] != 'NaN':
                        mean_volume_C1 = data['C1'].loc[data['C1']['Label'] == 'Mean', 'Volume (micron^3)']
                        sum_volume_C1 = data['C1'].iloc[:-4, 3].sum()
                elif data['C1'].shape[0] == 0:
                    mean_volume_C1 = 0
                    sum_volume_C1 = 0
                else:
                    mean_volume_C1 = data['C1']['Volume (micron^3)'][0]
                    sum_volume_C1 = data['C1']['Volume (micron^3)'][0]

                # Find mean volume, total volume and number of C2 elements
                C2_elements = data['C2']['Total C2'][0]
                if data['C2'].shape[0] > 1:
                    if data['C2']['Label'][0] != 'NaN':
                        mean_volume_C2 = data['C2'].loc[data['C2']['Label'] == 'Mean', 'Volume (micron^3)']
                        sum_volume_C2 = data['C2'].iloc[:-4, 3].sum()
                elif data['C2'].shape[0] == 0:
                    mean_volume_C2 = 0
                    sum_volume_C2 = 0
                else:
                    mean_volume_C2 = data['C2']['Volume (micron^3)'][0]
                    sum_volume_C2 = data['C2']['Volume (micron^3)'][0]

                # Find mean volume of contacts elements
                if data['Global_Coloc_(Co)'].shape[0] > 1:
                    contact = data['Global_Coloc_(Co)'].shape[0] - 4
                    mean_volume_contact = data['Global_Coloc_(Co)'].loc[
                        data['Global_Coloc_(Co)']['Label'] == 'Mean', 'Volume (micron^3)']
                elif data['Global_Coloc_(Co)'].shape[0] == 0:
                    contact = 0
                    mean_volume_contact = 0
                else:
                    contact = 1
                    mean_volume_contact = data['Global_Coloc_(Co)']['Volume (micron)']

                # Find C3 control channel values
                if 'Label.1' in data['Global_Coloc_(Co)'].columns:
                    if data['Global_Coloc_(Co)']['Count.1'].max() > 1:
                        controlled_coloc_number = data['Global_Coloc_(Co)']['Count.1'].max() - 4
                        controlled_coloc_number_per_um3 = controlled_coloc_number / volume_stack
                        controlled_coloc_mean_vol = list(data['Global_Coloc_(Co)'].loc[
                                                             data['Global_Coloc_(Co)'][
                                                                 'Label.1'] == 'Mean', 'Volume (micron^3).1'])
                        controlled_coloc_total_vol = data['Global_Coloc_(Co)'].loc[:data['Global_Coloc_(Co)'].loc[
                                                                                        data['Global_Coloc_(Co)'][
                                                                                            'Label.1'] == 'Mean'].index[
                                                                                        0] - 1,
                                                     'Volume (micron^3).1'].sum()
                    else:
                        controlled_coloc_number = (data['Global_Coloc_(Co)']['Volume (micron^3).1'] > 0).sum()
                        controlled_coloc_number_per_um3 = controlled_coloc_number / volume_stack
                        controlled_coloc_mean_vol = data['Global_Coloc_(Co)']['Volume (micron^3).1'][0]
                        controlled_coloc_total_vol = data['Global_Coloc_(Co)']['Volume (micron^3).1'][0]

                # Find sum, volume and intensity of C1 elements interacting with C2
                sum_int_vol_C1 = data['C1xCo'].iloc[:-1, 3].sum()
                count_C1_contacting_C2 = (data['C1xCo'].iloc[:-1, 3] > 0).sum()
                if 'Mean Gray value Globalcoloc area on C1' in data['C1xCo'].columns:
                    if 'Mean' in data['C1xCo']['Label Global_Coloc']:
                        intensity_of_C1_contacting_C2 = data['C1xCo'].loc[
                            data['C1xCo']['Label Global_Coloc'] == 'Mean', 'Mean Gray value Globalcoloc area on C1']
                    else:
                        intensity_of_C1_contacting_C2 = data['C1xCo']['Mean Gray value Globalcoloc area on C1'][0]

                # Find mean number of interaction of C1/C2 elements
                total_coloc = (data['C1xCo']['Nbs of Globalcoloc on C1']).sum()

                if not pd.isna(data['C1xCo']['C1 elements.1'][0]):  # Check if controlled coloc analysis was performed
                    # Find sum, volume and intensity of controlled C1 elements interacting with C2
                    count_controlled_C1_contacting_C2 = (data['C1xCo']['Total volume of Controlcoloc micron^3'] != 0).sum()
                    count_controlled_C2_contacting_C1 = (data['C2xCo']['Total volume of Controlcoloc micron^3'] != 0).sum()

                    C1_mean_interaction_volume_controlled = \
                        data['C1xCo'].loc[(data['C1xCo']['Total volume of Controlcoloc micron^3'] != 0)][
                            'Total volume of Controlcoloc micron^3'][:-1].mean()
                    C2_mean_interaction_volume_controlled = \
                        data['C2xCo'].loc[(data['C2xCo']['Total volume of Controlcoloc micron^3'] != 0)][
                            'Total volume of Controlcoloc micron^3'][:-1].mean()
                    C1_total_interaction_volume_controlled = \
                        data['C1xCo'].loc[(data['C1xCo']['Total volume of Controlcoloc micron^3'] != 0)][
                            'Total volume of Controlcoloc micron^3'][:-1].sum()
                    C2_total_interaction_volume_controlled = \
                        data['C2xCo'].loc[(data['C2xCo']['Total volume of Controlcoloc micron^3'] != 0)][
                            'Total volume of Controlcoloc micron^3'][:-1].sum()
                    mean_coloc_controlled_C1 = (data['C1xCo']['Nbs of Controlcoloc on C1']).sum()/C1_elements
                    mean_coloc_controlled_C2 = (data['C2xCo']['Nbs of Controlcoloc on C2']).sum()/C2_elements
                    controlled_coloc_in_total_coloc = (controlled_coloc_number/total_coloc) * 100

                if (data['C1xCo']['Label Control_Coloc'] == 'MGV_Controlcoloc:C1').sum() > 1:
                    mean_gray_value_coloc_controlled_C1 = list(data['C1xCo'].loc[
                                                                   data['C1xCo'][
                                                                       'Label Control_Coloc'] == 'Mean', 'Mean Gray value Controlcoloc area on C1'])
                elif (data['C1xCo']['Label Control_Coloc'] == 'MGV_Controlcoloc:C1').sum() == 1:
                    mean_gray_value_coloc_controlled_C1 = list(data['C1xCo'].loc[
                                                                   data['C1xCo'][
                                                                       'Label Control_Coloc'] == 'MGV_Controlcoloc:C1', 'Mean Gray value Controlcoloc area on C1'])

                if (data['C2xCo']['Label Control_Coloc'] == 'MGV_Controlcoloc:C2').sum() > 1:
                    mean_gray_value_coloc_controlled_C2 = list(data['C2xCo'].loc[
                                                                   data['C2xCo'][
                                                                       'Label Control_Coloc'] == 'Mean', 'Mean Gray value Controlcoloc area on C2'])
                elif (data['C2xCo']['Label Control_Coloc'] == 'MGV_Controlcoloc:C2').sum() == 1:
                    mean_gray_value_coloc_controlled_C2 = list(data['C2xCo'].loc[
                                                                   data['C2xCo'][
                                                                       'Label Control_Coloc'] == 'MGV_Controlcoloc:C2', 'Mean Gray value Controlcoloc area on C2'])

                # Find sum, volume and intensity of C2 elements interacting with C1
                sum_int_vol_C2 = data['C2xCo'].iloc[:-1, 3].sum()
                count_C2_contacting_C1 = (data['C2xCo'].iloc[:-1, 3] > 0).sum()
                if 'Mean Gray value Globalcoloc area on C2' in data['C2xCo'].columns:
                    if 'Mean' in data['C2xCo']['Label Global_Coloc']:
                        intensity_of_C2_contacting_C1 = data['C2xCo'].loc[
                            data['C2xCo']['Label Global_Coloc'] == 'Mean', 'Mean Gray value Globalcoloc area on C2']
                    else:
                        intensity_of_C2_contacting_C1 = data['C2xCo']['Mean Gray value Globalcoloc area on C2'][0]

                # Calculate the mean of interaction volume of C1 elements
                if C1_elements > 0:
                    mean_int_vol_C1 = sum_int_vol_C1 / C1_elements
                    percentage_C1_coloc_C2 = (count_C1_contacting_C2 / C1_elements) * 100
                else:
                    mean_int_vol_C1 = 0  # Isn't better to put NaN to not count this ?
                    percentage_C1_coloc_C2 = 0

                # Calculate the mean of interaction volume of C2 elements
                if C2_elements > 0:
                    mean_int_vol_C2 = sum_int_vol_C2 / C2_elements
                    percentage_C2_coloc_C1 = (count_C2_contacting_C1 / C2_elements) * 100
                else:
                    mean_int_vol_C2 = 0
                    percentage_C2_coloc_C1 = 0

                # Calculate number of C1/C2 elements per v.u
                nb_C1_vol = (C1_elements / volume_stack).values
                nb_C2_vol = (C2_elements / volume_stack).values

                data_to_plot[1] = pd.concat(
                    [data_to_plot[1], pd.Series(total_coloc / C1_elements, index=[condition[cond]['cond_name']],
                                                dtype='object')])  # Mean of coloc per C1 elements
                data_to_plot[2] = pd.concat(
                    [data_to_plot[2], pd.Series(total_coloc / C2_elements, index=[condition[cond]['cond_name']],
                                                dtype='object')])  # Mean of coloc per C2 elements
                data_to_plot[3] = pd.concat([data_to_plot[3], pd.Series(nb_C1_vol, index=[condition[cond]['cond_name']],
                                                                        dtype='object')])  # Nbs of C1 elements per v.u
                data_to_plot[4] = pd.concat([data_to_plot[4], pd.Series(nb_C2_vol, index=[condition[cond]['cond_name']],
                                                                        dtype='object')])  # Nbs of C2 elements per v.u
                if not np.isnan(np.sum(mean_volume_C1)):
                    data_to_plot[5] = pd.concat([data_to_plot[5],
                                                 pd.Series(mean_volume_C1.values, index=[condition[cond]['cond_name']],
                                                           dtype='object')])  # Mean volume of C1 elements
                    data_to_plot[6] = pd.concat([data_to_plot[6],
                                                 pd.Series(mean_volume_C2.values, index=[condition[cond]['cond_name']],
                                                           dtype='object')])  # Mean volume of C2 elements
                    data_to_plot[15] = pd.concat([data_to_plot[15],
                                                  pd.Series(sum_volume_C1, index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # Total volume of C1 elements
                    data_to_plot[16] = pd.concat([data_to_plot[16],
                                                  pd.Series(sum_volume_C2, index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # Total volume of C2 elements
                data_to_plot[7] = pd.concat([data_to_plot[7],
                                             pd.Series(mean_volume_contact.values, index=[condition[cond]['cond_name']],
                                                       dtype='object')])  # Mean of interaction volume
                data_to_plot[8] = pd.concat([data_to_plot[8],
                                             pd.Series(sum_int_vol_C1, index=[condition[cond]['cond_name']],
                                                       dtype='object')])  # Total interaction volume with C1 elements
                data_to_plot[9] = pd.concat([data_to_plot[9],
                                             pd.Series(sum_int_vol_C2, index=[condition[cond]['cond_name']],
                                                       dtype='object')])  # Total interaction volume with C2 elements
                data_to_plot[10] = pd.concat([data_to_plot[10], pd.Series((contact / volume_stack).values,
                                                                          index=[condition[cond]['cond_name']],
                                                                          dtype='object')])  # Total number of coloc per v.u
                data_to_plot[11] = pd.concat([data_to_plot[11],
                                              pd.Series(mean_int_vol_C1, index=[condition[cond]['cond_name']],
                                                        dtype='object')])  # Mean of interaction volume per C1 elements
                data_to_plot[12] = pd.concat([data_to_plot[12],
                                              pd.Series(mean_int_vol_C2, index=[condition[cond]['cond_name']],
                                                        dtype='object')])  # Mean of interaction volume per C2 elements
                data_to_plot[13] = pd.concat([data_to_plot[13],
                                              pd.Series(percentage_C1_coloc_C2, index=[condition[cond]['cond_name']],
                                                        dtype='object')])  # % C1 elements coloc with C2
                data_to_plot[14] = pd.concat([data_to_plot[14],
                                              pd.Series(percentage_C2_coloc_C1, index=[condition[cond]['cond_name']],
                                                        dtype='object')])  # % C2 elements coloc with C1

                data_to_plot[21] = pd.concat([data_to_plot[21], pd.Series(contact, index=[condition[cond]['cond_name']],
                                                                          dtype='object')])  # Total number of coloc
                data_to_plot[25] = pd.concat([data_to_plot[25],
                                              pd.Series(intensity_of_C1_contacting_C2,
                                                        index=[condition[cond]['cond_name']],
                                                        dtype='object')])  # MGV of C1 element contacting C2 element
                data_to_plot[26] = pd.concat([data_to_plot[26],
                                              pd.Series(intensity_of_C2_contacting_C1,
                                                        index=[condition[cond]['cond_name']],
                                                        dtype='object')])  # MGV of C2 element contacting C1 element

                if not 'Count.2' in data['C1'].columns:
                    data_to_plot[17] = pd.concat([data_to_plot[17],
                                                  pd.Series(nucleus_number_C1, index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # Nbs of nucleus coloc with C1
                    data_to_plot[18] = pd.concat([data_to_plot[18],
                                                  pd.Series(nucleus_number_C2, index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # Nbs of nucleus coloc with C2
                    data_to_plot[19] = pd.concat([data_to_plot[19],
                                                  pd.Series(nucleus_mean_area_C1, index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # C1 mean area coloc with nucleus
                    data_to_plot[20] = pd.concat([data_to_plot[20],
                                                  pd.Series(nucleus_mean_area_C2, index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # C2 mean area coloc with nucleus

                if 'Label.1' in data['Global_Coloc_(Co)'].columns:
                    data_to_plot[22] = pd.concat([data_to_plot[22],
                                                  pd.Series(controlled_coloc_number,
                                                            index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # Nbs of controlled coloc with C3
                    data_to_plot[23] = pd.concat([data_to_plot[23],
                                                  pd.Series(controlled_coloc_mean_vol,
                                                            index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # Mean volume of controlled coloc with C3
                    data_to_plot[24] = pd.concat([data_to_plot[24],
                                                  pd.Series(controlled_coloc_total_vol,
                                                            index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # Total volume of controlled coloc with C3
                    data_to_plot[37] = pd.concat([data_to_plot[37],
                                                  pd.Series(controlled_coloc_number_per_um3.values,
                                                            index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # Nbs of controlled coloc with C3 per v.u
                if not pd.isna(data['C1xCo']['C1 elements.1'][0]):
                    data_to_plot[27] = pd.concat([data_to_plot[27],
                                                  pd.Series((count_controlled_C1_contacting_C2 * 100 / C1_elements),
                                                            index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # % C1 elements coloc controlled with C2
                    data_to_plot[28] = pd.concat([data_to_plot[28],
                                                  pd.Series((count_controlled_C2_contacting_C1 * 100 / C2_elements),
                                                            index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # % C2 elements coloc controlled with C1
                    data_to_plot[29] = pd.concat([data_to_plot[29],
                                                  pd.Series(C1_mean_interaction_volume_controlled,
                                                            index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # C1 mean of interaction volume controlled
                    data_to_plot[30] = pd.concat([data_to_plot[30],
                                                  pd.Series(C2_mean_interaction_volume_controlled,
                                                            index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # C2 mean of interaction volume controlled
                    data_to_plot[31] = pd.concat([data_to_plot[31],
                                                  pd.Series(C1_total_interaction_volume_controlled,
                                                            index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # C1 total interaction volume controlled
                    data_to_plot[32] = pd.concat([data_to_plot[32],
                                                  pd.Series(C2_total_interaction_volume_controlled,
                                                            index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # C2 total interaction volume controlled
                    data_to_plot[33] = pd.concat([data_to_plot[33],
                                                  pd.Series(mean_coloc_controlled_C1,
                                                            index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # Mean of controlled coloc per C1 elements
                    data_to_plot[34] = pd.concat([data_to_plot[34],
                                                  pd.Series(mean_coloc_controlled_C2,
                                                            index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # Mean of controlled coloc per C2 elements
                    data_to_plot[35] = pd.concat([data_to_plot[35],
                                                  pd.Series(mean_gray_value_coloc_controlled_C1,
                                                            index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # MGV of C1 elements on controlled coloc site
                    data_to_plot[36] = pd.concat([data_to_plot[36],
                                                  pd.Series(mean_gray_value_coloc_controlled_C2,
                                                            index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # MGV of C2 elements on controlled coloc site
                    data_to_plot[38] = pd.concat([data_to_plot[38],
                                                  pd.Series(controlled_coloc_in_total_coloc,
                                                            index=[condition[cond]['cond_name']],
                                                            dtype='object')])  # % of controlled coloc in total coloc

    writer = pd.ExcelWriter(f'{save_path}/stat.xlsx', engine='xlsxwriter')
    data_writer = pd.ExcelWriter(f'{save_path}/data_excel.xlsx', engine='xlsxwriter')

    def estimation_plot(data_series, plot_title, outputpath):
        dataframe_estimation = pd.DataFrame(
            {'data': name_list, 'cond': data_series.index, 'value': data_series.values.astype(float)})
        data_estimation_plot = dabest.load(data=dataframe_estimation, x='cond', y='value', idx=cond_names)
        p = data_estimation_plot.mean_diff.plot(group_summaries='mean_sd',
                                                swarm_label=plot_title,
                                                custom_palette=colors,
                                                swarm_desat=.6,
                                                contrast_label='Mean difference',
                                                halfviolin_alpha=.5,
                                                halfviolin_desat=1)

        rawdata_axes = p.axes[0]
        effsize_axes = p.axes[1]
        rawdata_axes.set_xticklabels(cond_names)
        plt.savefig(outputpath, bbox_inches='tight')

        data_estimation_plot.mean_diff.statistical_tests.to_excel(writer, sheet_name=plot_title[:31])
        dataframe_estimation.to_excel(data_writer, sheet_name=plot_title[:31], index=name_list)

    estimation_plot(data_to_plot[1], plot_title='Mean of coloc per C1 elem',
                    outputpath=f'{C1xCo_save_path}/Mean of coloc per C1 elem.pdf')
    estimation_plot(data_to_plot[2], plot_title='Mean of coloc per C2 elem',
                    outputpath=f'{C2xCo_save_path}/Mean of coloc per C2 elem.pdf')
    estimation_plot(data_to_plot[3], plot_title='Nbs of C1 elem per v.u',
                    outputpath=f'{C1_save_path}/Nbs of C1 elem per v.u.pdf')
    estimation_plot(data_to_plot[4], plot_title='Nbs of C2 elem per v.u',
                    outputpath=f'{C2_save_path}/Nbs of C2 elem per v.u.pdf')
    if not np.isnan(np.sum(mean_volume_C1)):
        estimation_plot(data_to_plot[5], plot_title='Mean volume of C1 elem',
                        outputpath=f'{C1_save_path}/Mean volume of C1 elem.pdf')
        estimation_plot(data_to_plot[6], plot_title='Mean volume of C2 elem',
                        outputpath=f'{C2_save_path}/Mean volume of C2 elem.pdf')
    estimation_plot(data_to_plot[7], plot_title='Mean of i.v.',
                    outputpath=f'{Coloc_save_path}/Mean of i.v..pdf')
    estimation_plot(data_to_plot[8], plot_title='Total i.v. with C1 elem',
                    outputpath=f'{C1xCo_save_path}/Total i.v. with C1 elem.pdf')
    estimation_plot(data_to_plot[9], plot_title='Total i.v. with C2 elem',
                    outputpath=f'{C2xCo_save_path}/Total i.v. with C2 elem.pdf')
    estimation_plot(data_to_plot[10], plot_title='Total number of coloc per v.u',
                    outputpath=f'{Coloc_save_path}/Total number of coloc per v.u.pdf')
    estimation_plot(data_to_plot[11], plot_title='Mean of i.v. per C1 elem',
                    outputpath=f'{C1xCo_save_path}/Mean of i.v. per C1 elem.pdf')
    estimation_plot(data_to_plot[12], plot_title='Mean of i.v. per C2 elem',
                    outputpath=f'{C2xCo_save_path}/Mean of i.v. per C2 elem.pdf')
    estimation_plot(data_to_plot[13], plot_title='% C1 elem coloc with C2',
                    outputpath=f'{C1xCo_save_path}/% C1 elem coloc with C2.pdf')
    estimation_plot(data_to_plot[14], plot_title='% C2 elem coloc with C1',
                    outputpath=f'{C2xCo_save_path}/% C2 elem coloc with C1.pdf')
    if not np.isnan(np.sum(mean_volume_C1)):
        estimation_plot(data_to_plot[15], plot_title='Total volume of C1 elem',
                        outputpath=f'{C1_save_path}/Total volume of C1 elem.pdf')
        estimation_plot(data_to_plot[16], plot_title='Total volume of C2 elem',
                        outputpath=f'{C2_save_path}/Total volume of C2 elem.pdf')
    estimation_plot(data_to_plot[21], plot_title='Total number of coloc',
                    outputpath=f'{Coloc_save_path}/Total number of coloc.pdf')
    estimation_plot(data_to_plot[25], plot_title='MGV of C1 elem on coloc site',
                    outputpath=f'{C1xCo_save_path}/MGV of C1 elem on coloc site.pdf')
    estimation_plot(data_to_plot[26], plot_title='MGV of C2 elem on coloc site',
                    outputpath=f'{C2xCo_save_path}/MGV of C2 elem on coloc site.pdf')

    if not 'Count.2' in data['C1'].columns:
        estimation_plot(data_to_plot[17], plot_title='Number of nucleus in C1',
                        outputpath=f'{C1_save_path}/Number of nucleus in C1.pdf')
        estimation_plot(data_to_plot[18], plot_title='Number of nucleus in C2',
                        outputpath=f'{C2_save_path}/Number of nucleus in C2.pdf')
        estimation_plot(data_to_plot[19], plot_title='C1 soma mean area',
                        outputpath=f'{C1_save_path}/C1 soma mean area.pdf')
        estimation_plot(data_to_plot[20], plot_title='C2 soma mean area',
                        outputpath=f'{C2_save_path}/C2 soma mean area.pdf')

    if 'Label.1' in data['Global_Coloc_(Co)'].columns:
        estimation_plot(data_to_plot[22], plot_title='Number of ctld coloc C3',
                        outputpath=f'{Coloc_save_path}/Number of ctld coloc C3.pdf')
        estimation_plot(data_to_plot[23], plot_title='Mean volume of ctld coloc C3',
                        outputpath=f'{Coloc_save_path}/Mean volume of ctld coloc C3.pdf')
        estimation_plot(data_to_plot[24], plot_title='Total volume of ctld coloc C3',
                        outputpath=f'{Coloc_save_path}/Total volume of ctld coloc C3.pdf')
        estimation_plot(data_to_plot[37], plot_title='Total number of ctld coloc per v.u',
                        outputpath=f'{Coloc_save_path}/Total number of ctld coloc per v.u.pdf')

    if not pd.isna(data['C1xCo']['C1 elements.1'][0]):
        estimation_plot(data_to_plot[27], plot_title='% C1 elem coloc ctld with C2',
                        outputpath=f'{C1xCo_save_path}/% C1 elem coloc ctld with C2 .pdf')
        estimation_plot(data_to_plot[28], plot_title='% C2 elem coloc ctld with C1',
                        outputpath=f'{C2xCo_save_path}/% C2 elem coloc ctld with C1.pdf')
        estimation_plot(data_to_plot[29], plot_title='Mean of ctld i.v. per C1 elem',
                        outputpath=f'{C1xCo_save_path}/Mean of ctld i.v. per C1 elem.pdf')
        estimation_plot(data_to_plot[30], plot_title='Mean of ctld i.v. per C2 elem',
                        outputpath=f'{C2xCo_save_path}/Mean of ctld i.v. per C2 elem.pdf')
        estimation_plot(data_to_plot[31], plot_title='Total ctld i.v. with C1 elem',
                        outputpath=f'{C1xCo_save_path}/Total ctld i.v. with C1 elem.pdf')
        estimation_plot(data_to_plot[32], plot_title='Total ctld i.v. with C2 elem',
                        outputpath=f'{C2xCo_save_path}/Total ctld i.v. with C2 elem.pdf')
        estimation_plot(data_to_plot[33], plot_title='Mean of ctld coloc per C1 elem',
                        outputpath=f'{C1xCo_save_path}/Mean of ctld coloc per C1 elem.pdf')
        estimation_plot(data_to_plot[34], plot_title='Mean of ctld coloc per C2 elem',
                        outputpath=f'{C2xCo_save_path}/Mean of ctld coloc per C2 elem.pdf')
        estimation_plot(data_to_plot[35], plot_title='MGV of C1 elem on ctld coloc site',
                        outputpath=f'{C1xCo_save_path}/MGV of C1 elem on ctld coloc site.pdf')
        estimation_plot(data_to_plot[36], plot_title='MGV of C2 elem on ctld coloc site',
                        outputpath=f'{C2xCo_save_path}/MGV of C2 elem on ctld coloc site.pdf')
        estimation_plot(data_to_plot[38], plot_title='% of ctld coloc in total coloc',
                        outputpath=f'{Coloc_save_path}/% of ctld coloc in total coloc.pdf')
    writer.close()
    data_writer.close()

class App(ctk.CTk):
    def __init__(self):
        super().__init__()

        global list_entry, list_color_entry
        list_entry = []
        list_color_entry = []

        # Create the main window
        self.title("Bioloc3D")
        self.after(201, lambda: self.iconbitmap('./ressource/logo.ico'))
        self.grid_columnconfigure(2, weight=1)

        # Create a sidebar at the left of the main window
        self.sidebar_frame = ctk.CTkFrame(self, width=140, corner_radius=0)
        self.sidebar_frame.grid(row=0, column=0, rowspan=20, sticky="ns")
        self.sidebar_frame.grid_rowconfigure(2, weight=1)
        # Add a title to the sidebar
        self.logo_label = ctk.CTkLabel(self.sidebar_frame, text="Bioloc3D", font=ctk.CTkFont(size=20, weight="bold"))
        self.logo_label.grid(row=0, column=0, padx=40, pady=(10, 0), sticky="nw")

        # Add logo image on the sidebar
        self.image = ctk.CTkLabel(self.sidebar_frame,
                                  image=ctk.CTkImage(Image.open(logo_png),
                                                     size=(125, 125)),
                                  text=None)
        self.image.grid(row=1, column=0, padx=20, pady=(5, 0), sticky="n")

        # Create a progress bar
        self.progress_bar = ctk.CTkProgressBar(self.sidebar_frame, width=125, mode="indeterminate")
        self.progress_bar.grid(row=4, column=0, padx=20, pady=(10, 10), sticky="s")

        # Add a button to the sidebar permitting modification of the appearance mode
        self.switch_variable = ctk.StringVar(value='on')
        self.button = ctk.CTkSwitch(self.sidebar_frame, text="Dark mode", command=self.switch_appearance_mode_event,
                                    variable=self.switch_variable, onvalue="on", offvalue="off")
        self.button.grid(row=5, column=0, padx=20, pady=(10, 10), sticky="s")

        # Create button and label to open and save files
        self.load_file_button = ctk.CTkButton(self, text="Load file", command=self.load_file)
        self.load_file_button.grid(row=0, column=1, padx=5, pady=(5, 0), sticky="nw")

        self.load_file_label = ctk.CTkLabel(self, text="No file selected")
        self.load_file_label.grid(row=1, column=1, padx=5, pady=(0, 5), columnspan=2, sticky="nw")

        self.save_file_button = ctk.CTkButton(self, text="Save file", command=self.save_file)
        self.save_file_button.grid(row=2, column=1, padx=5, pady=(0, 0), sticky="nw")

        self.save_file_label = ctk.CTkLabel(self, text="No file selected")
        self.save_file_label.grid(row=3, column=1, padx=10, pady=(0, 0), columnspan=2, sticky="nw")

        # Create a label and combobox to select the number of condition
        self.number_condition_box = ctk.CTkComboBox(self, values=["2", "3", "4", "5", "6", "7", "8", "9", "10"],
                                                    command=self.combobox_callback, state="readonly", width=155)
        self.number_condition_box.set('Number of condition')
        self.number_condition_box.grid(row=4, column=1, padx=5, pady=(10, 10), sticky="w")

        # Create a button to launch the analysis
        self.launch_button = ctk.CTkButton(self.sidebar_frame, text="Launch analysis", command=self.launch_analysis)
        self.launch_button.grid(row=3, column=0, padx=5, pady=(10, 10), sticky="s")

        self.protocol("WM_DELETE_WINDOW", self.quit)

    def switch_appearance_mode_event(self):
        if self.switch_variable.get() == 'on':
            ctk.set_appearance_mode('Dark')
        else:
            ctk.set_appearance_mode('Light')

    def load_file(self):
        self.load_filename = filedialog.askdirectory()
        self.load_file_label.configure(text=self.load_filename)

    def save_file(self):
        self.save_filename = filedialog.askdirectory()
        self.save_file_label.configure(text=self.save_filename)

    def combobox_callback(self, choice):
        self.number_of_condition = int(choice)
        self.entry_names = []
        self.colors = []

        self.colors_label = ctk.CTkLabel(self, text="Colors \n(seaborn or hexadecimal)")
        self.colors_label.grid(row=4, column=2, padx=6, pady=(0, 0), sticky="w")
        for x in list_entry:
            x.destroy()
        for x in list_color_entry:
            x.destroy()
        for i in range(1, self.number_of_condition + 1):
            self.entry_loop = ctk.CTkEntry(self, placeholder_text="Condition {}".format(i))
            self.entry_loop.grid(row=4 + i, column=1, padx=5, pady=(5, 5), sticky="w")
            list_entry.append(self.entry_loop)
            self.entry_names.append(self.entry_loop)

            self.color_entry_loop = ctk.CTkEntry(self, placeholder_text="Color {}".format(i))
            self.color_entry_loop.grid(row=4 + i, column=2, padx=5, pady=(5, 5), sticky="w")
            list_color_entry.append(self.color_entry_loop)
            self.colors.append(self.color_entry_loop)

    def launch_analysis(self):
        self.cond_name = [x.get() for x in self.entry_names]
        self.color_palette = [x.get() for x in self.colors]
        self.progress_bar.start()

        # Lancer le script de quantification dans un thread séparé
        analysis_thread = threading.Thread(target=self.run_quantification_script)
        analysis_thread.start()

    def run_quantification_script(self):
        try:
            quantification_script(
                pathway_to_raw_data=self.load_filename,
                nb_cond=self.number_of_condition,
                cond_names=self.cond_name,
                save_path=self.save_filename,
                colors=self.color_palette if self.color_palette[0] != "" else None
            )
            # Planifier l'affichage du message dans le thread principal
            self.after(0, self.show_completion_message)
        except Exception as e:
            # Planifier l'affichage du message d'erreur dans le thread principal
            self.after(0, lambda: self.show_error_message(str(e)))
        finally:
            # Planifier l'arrêt de la barre de progression dans le thread principal
            self.after(0, self.progress_bar.stop)

    def show_completion_message(self):
        CTkMessagebox(title="Analysis finished",
                      message="The analysis is finished, you can find the results in the folder you selected")

    def show_error_message(self, error_message):
        CTkMessagebox(title="Error", message=error_message)

if __name__ == "__main__":
    app = App()
    app.mainloop()
