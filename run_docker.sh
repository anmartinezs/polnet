# Read the arguments for the config yaml file and outdir
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --config)
            config_file="$2"
            shift # past argument
            shift # past value
            ;;
        --out_dir)
            out_dir="$2"
            shift # past argument
            shift # past value
            ;;
        *)    # unknown option
    esac
done

# Check if the config file argument is provided
if [[ -z $config_file ]]; then
    echo "Error: Config file argument (--config) is missing."
    exit 1
fi

# Check if the outdir argument is provided
if [[ -z $out_dir ]]; then
    echo "Error: Outdir argument (--out_dir) is missing."
    exit 1
fi

# prints ids
echo "myUID: $(id -u)"
echo "myGID: $(id -g)"

# Run the docker container with the provided config file and outdir
docker run -e PYTHONUNBUFFERED=1 -v "$config_file":/app/generation_config.yaml -v "$out_dir":/app/outdir --user $(id -u):$(id -g) polnet_docker
