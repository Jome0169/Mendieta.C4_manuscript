import argparse
import pybedtools
import csv

def modify_interval(interval):
    """
    Function to modify a pybedtools.Interval object.
    Keeps the first three columns and adds a fourth column with chrom_start_stop format.
    """
    # Create the new name for the fourth column
    new_name = f"{interval.chrom}_{interval.start}_{interval.end}"
    # Modify the interval to only have the first three columns and the new fourth column
    interval.name = new_name
    return interval



def reformat(x):
    for line in open(x.fn):
        fields = line.strip().split('\t')
        yield pybedtools.create_interval_from_list([
            fields[0], fields[1], fields[2]])



def integrate_sites(regions_bed, tn5_sites, output_file):
    # Load the regions from the BED file using pybedtools
    regions = pybedtools.BedTool(regions_bed)



    updated_region_names = regions.each(modify_interval)

    
    # Load the Tn5 integration sites directly with pybedtools, without using pandas
    tn5_sites_bed = pybedtools.BedTool(tn5_sites)
    
    # Intersect the two sets, keeping information about the Tn5 integrations
    intersections = updated_region_names.intersect(tn5_sites_bed, wa=True, wb=True)

    #$print(intersections)
    


    region_BC_dict = {}
    test_group = intersections.groupby(g="4,9", c="5", o="count").saveas()

    for line in open(test_group.fn):
        line = line.strip().split('\t')
        chrom, barcode = line[0], line[1]
    
        # Initialize the nested dictionary for 'chrom' if it doesn't exist
        if chrom not in region_BC_dict:
            region_BC_dict[chrom] = {}
    
        # Increment the count for the barcode, initializing it if it doesn't exist
        region_BC_dict[chrom][barcode] = region_BC_dict[chrom].get(barcode, 1)
        

    # Open the file for writing
    with open(output_file, "w", newline='') as file:
        # Create a CSV writer object for writing to a TSV file
        writer = csv.writer(file, delimiter='\t')
        
        # Write the header row with row_number as the first column
        writer.writerow(["", "i", "j", "x"])
        
        # Initialize a counter for the row number
        row_number = 1
        
        # Iterate over each item in the dictionary
        for region, bc_vals_list in region_BC_dict.items():
            for bc_val, value in bc_vals_list.items():
                # Write the row number, region, BC_val, and its value to the file, with row_number first
                writer.writerow([row_number, region, bc_val, value])
                # Increment the row number for the next entry
                row_number += 1


def main():
    parser = argparse.ArgumentParser(description='Integrate Tn5 sites with genomic regions.')
    parser.add_argument('-regions_bed', required=True, help='Path to BED file containing regions')
    parser.add_argument('-tn5_sites', required=True, help='Path to file containing Tn5 integration sites')
    parser.add_argument('-output_file', required=True, help='Path to save the output')
    args = parser.parse_args()
    integrate_sites(args.regions_bed, args.tn5_sites, args.output_file)

if __name__ == '__main__':
    main()
