import psutil
import time
import datetime
import sys


def monitor_processes(pid, keyword, log_file, interval = 1):
    with open(log_file, "a") as file:
        while True:
            try:
                process = psutil.Process(pid)
            except psutil.NoSuchProcess:
                return
            matching_processes = []
            
            # Find all processes matching the keyword
            for process in psutil.process_iter(attrs=['pid', 'name', 'cmdline']):
                try:
                    name = process.info['name']
                    cmdline = ' '.join(process.info['cmdline']) if process.info['cmdline'] else ''
                    if keyword in name or keyword in cmdline:
                        matching_processes.append(process)
                        # Take an initial CPU usage snapshot
                        process.cpu_percent(interval=None)  # Start tracking CPU usage
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    continue
                
            if len(matching_processes) == 1:
                if 'monitor_pid.py' in cmdline:
                    break
                
            # Wait for the sampling interval
            psutil.time.sleep(interval)

            # Aggregate CPU and memory stats
            total_cpu = 0.0
            total_rss = 0
            # total_vms = 0
                    
            for process in matching_processes:
                try:
                    total_cpu += process.cpu_percent(interval=None)  # Get CPU percent
                    mem_info = process.memory_info()  # Get memory info
                    total_rss += mem_info.rss / (1024 * 1024) 
                    # total_vms += mem_info.vms
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    continue
            # Get current timestamp
            timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            # Log the information to the file
            file.write(f"{timestamp} - Num Process: {len(matching_processes)} - RAM Usage: {total_rss:.2f} MB  - CPU Usage: {total_cpu:.2f}%\n")
            file.flush()  # Flush to ensure data is written immediately

            time.sleep(log_interval) 
    
def log_ram_cpu_usage(pid, log_file):
    try:
        # Get the process object by PID
        process = psutil.Process(pid)
        with open(log_file, "a") as file:
            while True:
                # Get current memory usage in MB
                memory_info = process.memory_info()
                ram_usage_mb = memory_info.rss / (1024 * 1024)  # Convert bytes to MB

                # Get current CPU usage as a percentage
                cpu_usage_percent = process.cpu_percent(interval=1)  # 1 second sample interval

                # Get current timestamp
                timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

                # Log the information to the file
                file.write(f"{timestamp} - PID {pid} - RAM Usage: {ram_usage_mb:.2f} MB - CPU Usage: {cpu_usage_percent:.2f}%\n")
                file.flush()  # Flush to ensure data is written immediately

                # Wait for the specified interval before logging again
                time.sleep(log_interval - 1)  # Adjust sleep since cpu_percent already waited 1 second
    except psutil.NoSuchProcess:
        print(f"Process with PID {pid} does not exist.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    pid = int(sys.argv[1])  # Example PID, replace with your actual PID
    keyword = sys.argv[2] #'MaxQuant'
    log_file = sys.argv[3]  # Log file path
    log_interval = 10  # Time interval in seconds between logs

    monitor_processes(pid, keyword, log_file)
    # log_ram_cpu_usage(pid, log_file)
